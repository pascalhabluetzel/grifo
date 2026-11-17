"""
cluster.py — reference-free clustering of amplicon reads into consensus
sequences, plus read-to-cluster assignment.

cluster()      groups reads into clusters and returns one consensus
               sequence per cluster
assign_reads() assigns every read to its best-matching cluster and counts
               reads per cluster

Works with pandas 1.x and 2.x, Python >= 3.8. Requires the external tools
minimap2 and spoa on your PATH (or pass their paths with minimap2= / spoa=).

Usage
-----
    from cluster import cluster, assign_reads
    centers = cluster('barcode01_concatenated.csv',
                      'barcode01_clusters.csv',
                      n_iter=5)
    assignments, counts = assign_reads('barcode01_concatenated.csv',
                                       'barcode01_clusters.csv')

Input : csv with columns [id, sequence] (other columns are ignored), fasta,
        or fastq.
Output: csv with columns [id, sequence, split]; ids are cluster0, cluster1, ...

License
-------
The clustering algorithm is a port of code that is
Copyright (C) 2020 Bilgenur Baloglu, Centre for Biodiversity Genomics,
University of Guelph, distributed under the GNU Lesser General Public
License version 3. As a derivative work, this file is distributed under the
same license (https://www.gnu.org/licenses/lgpl-3.0.html).
"""

import os
import re
import shutil
import logging
import gzip
import tempfile
import warnings
import subprocess

import numpy as np
import pandas as pd
import sklearn.cluster

logger = logging.getLogger('cluster')

# default clustering parameters
DEFAULTS = dict(
    n_iter=10,                           # -iter   clst_N_iter
    csize=20,                            # -cs     clst_csize
    N=2000,                              # -N      clst_N
    th_s=4,                              # -ts     clst_th_s
    th_m=0.9,                            # -tm     clst_th_m
    pw_config='-k15 -w10 -p 0.9 -D',     # -pw_config clst_pw_config
    msa_config='-n -15 -g -10 -l 0 -r 0' # msa_config (spoa)
)


#############################################################################
# Main entry point
#############################################################################

def cluster(input_file, output_file=None, n_iter=DEFAULTS['n_iter'],
            csize=DEFAULTS['csize'], N=DEFAULTS['N'], th_s=DEFAULTS['th_s'],
            th_m=DEFAULTS['th_m'], pw_config=DEFAULTS['pw_config'],
            msa_config=DEFAULTS['msa_config'], init_centers=None,
            minimap2='minimap2', spoa='spoa', threads=3, seed=None,
            workspace=None, log_file=None):
    '''
    Cluster sequences into cluster-center consensus sequences.

    input_file   = csv with columns [id, sequence], or a fasta file
    output_file  = csv to write [id, sequence, split]; None to skip writing
    n_iter       = clustering iterations
    csize        = sequences per cluster used to build the spoa consensus
    N            = subsample size used for pairwise distances
    th_s         = number of k-means partitions in the sweep step
    th_m         = OPTICS eps threshold for merging similar clusters
    pw_config    = minimap2 options for read-vs-center alignment
    msa_config   = spoa options for consensus building
    init_centers = optional csv/fasta of starting cluster centers
    minimap2, spoa = executable names or paths
    threads      = threads per minimap2 call, and max parallel spoa jobs
    seed         = integer for reproducible runs (otherwise results vary
                   slightly between runs, as reads are sampled at random)
    workspace    = folder for temporary files (default: temp folder next to
                   output_file, removed afterwards)
    log_file     = optional log file path
    returns pandas DataFrame of cluster centers [id, sequence, split]
    '''
    _check_tool(minimap2)
    _check_tool(spoa)
    if seed is not None:
        np.random.seed(seed)
    handler = None
    if log_file:
        handler = logging.FileHandler(log_file, mode='w')
        handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s'))
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)

    out_dir = os.path.dirname(os.path.abspath(output_file)) if output_file else os.getcwd()
    own_workspace = workspace is None
    if own_workspace:
        workspace = tempfile.mkdtemp(prefix='cluster_', dir=out_dir)
    else:
        os.makedirs(workspace, exist_ok=True)

    tools = dict(minimap2=minimap2, spoa=spoa, threads=int(threads),
                 seed=seed, workspace=workspace)
    try:
        logger.info('Loading %s', input_file)
        df = _load_sequences(input_file)
        df['id'] = ['read' + str(i) for i in range(len(df))]
        # keep only reads with both primers if that info exists
        if 'fwd_primer_seq' in df.columns and 'rev_primer_seq' in df.columns:
            df = df[df['fwd_primer_seq'].notna() & df['rev_primer_seq'].notna()]
        df = df[['id', 'sequence']].reset_index(drop=True)
        logger.info('%d sequences to cluster', len(df))

        df_d = []
        if init_centers:
            df_d = _load_sequences(init_centers)[['id', 'sequence']]

        data = perform_cluster(df, df_d, max_iter=int(n_iter), csize=int(csize),
                               N=int(N), th_s=int(th_s), th_m=float(th_m),
                               pw_config=pw_config, msa_config=msa_config,
                               tools=tools)
        data = data[['id', 'sequence', 'split']].reset_index(drop=True)
        if output_file:
            data.to_csv(output_file, index=False)
            logger.info('Wrote %d clusters to %s', len(data), output_file)
        return data
    finally:
        if own_workspace:
            shutil.rmtree(workspace, ignore_errors=True)
        if handler is not None:
            logger.removeHandler(handler)
            handler.close()


def assign_reads(reads_file, clusters_file, assignments_file='auto',
                  counts_file='auto', min_identity=0.85, min_coverage=0.8,
                  minimap2='minimap2', threads=3, preset='map-ont', workspace=None):
    '''
    Assign every read to its best-matching cluster center and count reads per
    cluster. cluster() returns centers only; this adds read membership.

    reads_file    = the reads that were clustered: csv [id, sequence], fasta or
                    fastq (optionally .gz)
    clusters_file = output of cluster(): csv [id, sequence, ...] or fasta
    assignments_file = csv for the per-read table; 'auto' derives the name from
                    clusters_file (X_clusters.csv -> X_read_assignments.csv);
                    None to skip writing
    counts_file   = csv for reads per cluster; 'auto' -> X_cluster_counts.csv;
                    None to skip writing
    min_identity  = minimum alignment identity (matching bases / alignment
                    length) for a read to be assigned. Raw nanopore reads vs
                    a consensus typically score 0.85-0.97 depending on chemistry
    min_coverage  = minimum fraction of the read covered by the alignment
    preset        = minimap2 preset (-x)

    Each read goes to the cluster with the highest alignment score (AS).
    'ambiguous' = True when another cluster scored exactly as well; the read
    is still assigned (to the first such cluster) but you may want to
    exclude or split these for closely related clusters.

    returns (assignments, counts) as pandas DataFrames
      assignments: read_id, cluster, identity, coverage, AS, ambiguous, status
                   status is 'assigned', 'below_threshold' or 'no_hit'
      counts:      cluster, n_reads, n_ambiguous, plus an 'unassigned' row
    '''
    _check_tool(minimap2)
    reads = _load_sequences(reads_file)
    # keep the first word of the header as the read id (drops '>' and the
    # 'runid=... ch=...' fields that nanopore headers carry)
    reads['read_id'] = [str(i).lstrip('>').split()[0] if str(i).strip() else 'read%d' % k
                        for k, i in enumerate(reads['id'])]
    clusters = _load_sequences(clusters_file)[['id', 'sequence']]
    logger.info('assign_reads: %d reads vs %d clusters', len(reads), len(clusters))

    own_workspace = workspace is None
    if own_workspace:
        base_dir = os.path.dirname(os.path.abspath(clusters_file))
        workspace = tempfile.mkdtemp(prefix='assign_reads_', dir=base_dir)
    tools = dict(minimap2=minimap2, threads=int(threads), workspace=workspace)
    try:
        aln = run_minimap2(reads[['read_id', 'sequence']].rename(columns={'read_id': 'id'}),
                           clusters, config='-x ' + preset, cigar=True, tools=tools)
    finally:
        if own_workspace:
            shutil.rmtree(workspace, ignore_errors=True)

    if len(aln) > 0:
        aln['identity'] = aln['match'] / aln['tot']
        aln['coverage'] = (aln['q_end'] - aln['q_start']) / aln['q_len']
        # best alignment per read; detect ties with a different cluster
        aln = aln.sort_values(['query_id', 'AS', 'identity'], ascending=[True, False, False])
        best = aln.drop_duplicates('query_id', keep='first').set_index('query_id')
        top_as = best['AS'].reindex(aln['query_id']).values
        ties = aln[(aln['AS'].values == top_as)].groupby('query_id')['database_id'].nunique()
        best['ambiguous'] = ties.reindex(best.index).fillna(1).values > 1
    else:
        best = pd.DataFrame(columns=['database_id', 'identity', 'coverage', 'AS', 'ambiguous'])

    out = pd.DataFrame({'read_id': reads['read_id'].values})
    b = best.reindex(out['read_id'])
    out['cluster'] = b['database_id'].values
    out['identity'] = b['identity'].values.astype(float).round(4)
    out['coverage'] = b['coverage'].values.astype(float).round(4)
    out['AS'] = b['AS'].values
    out['ambiguous'] = b['ambiguous'].fillna(False).values.astype(bool)
    hit = out['cluster'].notna().values
    ok = hit & (out['identity'].values >= min_identity) & (out['coverage'].values >= min_coverage)
    out['status'] = np.where(ok, 'assigned', np.where(hit, 'below_threshold', 'no_hit'))
    out.loc[~ok, 'cluster'] = None
    out.loc[~ok, 'ambiguous'] = False

    assigned = out[out['status'] == 'assigned']
    counts = pd.DataFrame({'cluster': clusters['id'].values})
    counts['n_reads'] = counts['cluster'].map(assigned['cluster'].value_counts()).fillna(0).astype(int)
    counts['n_ambiguous'] = counts['cluster'].map(
        assigned[assigned['ambiguous']]['cluster'].value_counts()).fillna(0).astype(int)
    counts = counts.sort_values('n_reads', ascending=False).reset_index(drop=True)
    counts = pd.concat([counts, pd.DataFrame({'cluster': ['unassigned'],
                                              'n_reads': [int((~ok).sum())],
                                              'n_ambiguous': [0]})], ignore_index=True)

    stem = re.sub(r'(_clusters)?\.(csv|fa|fasta)(\.gz)?$', '', clusters_file)
    if assignments_file == 'auto':
        assignments_file = stem + '_read_assignments.csv'
    if counts_file == 'auto':
        counts_file = stem + '_cluster_counts.csv'
    if assignments_file:
        out.to_csv(assignments_file, index=False)
    if counts_file:
        counts.to_csv(counts_file, index=False)
    logger.info('assign_reads: %d assigned (%d ambiguous), %d below threshold, %d no hit',
                ok.sum(), out['ambiguous'].sum(), (out['status'] == 'below_threshold').sum(),
                (out['status'] == 'no_hit').sum())
    return out, counts


#############################################################################
# Clustering logic
#############################################################################

def perform_cluster(df_q, df_d, max_iter, csize, N, th_s, th_m, pw_config,
                    msa_config, tools):
    '''Refine cluster centers: initial sampling, then split / sweep / merge.'''
    if len(df_d) > 0:
        df_c = df_d.copy()
    else:
        df_c = cluster_sample(df_q, N=5, th=0.8, rsize=1, config=pw_config, tools=tools)
    df_c['split'] = True
    for k in range(max_iter):
        logger.info('perform_cluster: iter = %d/%d', k, max_iter)
        df_c = cluster_split(df_q, df_c, N, csize, pw_config, msa_config, tools)
        df_c = cluster_sweep(df_q, df_c, th_s, N, csize, pw_config, msa_config, tools)
        df_c = cluster_merge(df_q, df_c, th_m, 100, csize, pw_config, msa_config, tools)
    # final refinement
    df_c['split'] = True
    df_c = cluster_split(df_q, df_c, 200, csize, pw_config, msa_config, tools)
    df_c = cluster_merge(df_q, df_c, th_m, 100, csize, pw_config, msa_config, tools)
    logger.info('perform_cluster: complete, %d clusters', len(df_c))
    return df_c


def cluster_sample(df_q, N=5, th=0.8, rsize=1, config='', tools=None):
    '''Pick up to N rough, mutually dissimilar starting centers.'''
    qry = np.unique(df_q['id'].astype(str))
    if len(qry) == 0:
        logger.info('cluster_sample: nothing to cluster')
        return pd.DataFrame([], columns=['id', 'sequence'])
    # seed with a random read (the original code meant to do this, but its
    # DataFrame.append result was discarded, so the first round was wasted)
    ridx = np.random.permutation(len(qry))
    df_c = df_q[df_q['id'].isin(qry[ridx[:rsize]])][['id', 'sequence']]
    db = df_c
    while N > len(df_c):
        logger.info('cluster_sample: qlen=%d progress=%d/%d', len(qry), len(df_c), N)
        seq = df_q[df_q['id'].isin(qry)]
        out = run_minimap2(seq, db, config=config, cigar=True, tools=tools)
        out = out[out['match_score'] > th]
        qry = set(qry) - set(out['query_id'].astype(str)) - set(df_c['id'].astype(str))
        if len(qry) == 0:
            break
        qry = np.array(sorted(qry))
        ridx = np.random.permutation(len(qry))
        db = df_q[df_q['id'].isin(qry[ridx[:rsize]])][['id', 'sequence']]
        df_c = _concat([df_c, db])
    df_c = df_c.reset_index(drop=True)
    df_c['id'] = ['cluster' + str(i) for i in range(len(df_c))]
    return df_c


def cluster_eval(df_q, df_c, pw_config, tools, cigar=True):
    '''Align reads to centers; drop centers that attract no reads.'''
    df_align = run_minimap2(df_q, df_c, config=pw_config + ' --for-only',
                            cigar=cigar, tools=tools)
    cols = np.unique(df_align['database_id'])
    df_c = df_c[df_c['id'].isin(cols)]
    logger.info('cluster_eval: number of clusters = %d', len(df_c))
    return df_c, df_align


def cluster_split(df_q, df_c, N, csize, pw_config, msa_config, tools):
    '''Re-cluster the reads of each center marked split=True.'''
    df_c, df_align = cluster_eval(df_q, df_c, pw_config, tools, cigar=False)
    split = df_c['split'].astype(bool).values
    squeue = df_c['id'].values[split]
    df_c = df_c[~split]
    if len(df_align) > 0:
        df_align['weight'] = cluster_weight(df_align, metric='s1')
    for i, qout in enumerate(squeue):
        logger.info('cluster_split: splitting on cid=%s %d/%d', qout, i, len(squeue))
        qlist = cluster_subsample(df_align, [qout], N=N, mode='semi', metric='s1')
        qlist = df_q[df_q['id'].isin(qlist)]
        cout = cluster_compute(qlist, csize, pw_config=pw_config,
                               msa_config=msa_config, tools=tools)
        cout['split'] = len(cout) > 1
        df_c = _concat([df_c, cout[['id', 'sequence', 'split']]])
    return _rename_clusters(df_c)


def cluster_sweep(df_q, df_c, th_s, N, csize, pw_config, msa_config, tools):
    '''Find new centers among reads that align poorly to existing centers.'''
    df_c, df_align = cluster_eval(df_q, df_c, pw_config, tools, cigar=False)
    covered = set()
    if len(df_align) > 0:
        # normalise chaining score s1 by the best s1 for that center
        df_a = get_best(df_align, ['database_id'], metric='s1').rename(columns={'s1': 'b1'})
        df_align = df_align.merge(df_a[['database_id', 'b1']], on='database_id', how='left')
        df_align['m1'] = df_align['s1'] / df_align['b1']
        df_align = get_best(df_align, ['query_id'], metric='m1')
        df_align = df_align.rename(columns={'query_id': 'id'})
        if len(df_align) >= th_s:
            labels = cluster_kmeans_1d(df_align['m1'].values.astype(float), n_clusters=th_s,
                                       seed=tools['seed'])
            covered = set(df_align['id'].astype(str).values[labels != 0])  # 0 = worst partition
        else:
            covered = set(df_align['id'].astype(str))
    uncover = set(df_q['id'].astype(str)) - covered
    uncover = df_q[df_q['id'].isin(list(uncover))]
    logger.info('cluster_sweep: uncovered %d/%d', len(uncover), len(df_q))
    ridx = np.random.permutation(len(uncover))
    x = uncover.iloc[ridx[:N]]
    cout = cluster_compute(x, csize, pw_config=pw_config, msa_config=msa_config, tools=tools)
    cout['split'] = True
    df_c = _concat([df_c, cout[['id', 'sequence', 'split']]])
    return _rename_clusters(df_c)


def cluster_merge(df_q, df_c, th_m, N, csize, pw_config, msa_config, tools):
    '''Merge centers whose read-alignment profiles are similar.'''
    df_c, df_align = cluster_eval(df_q, df_c, pw_config, tools, cigar=True)
    if len(df_align) == 0:
        return _rename_clusters(df_c)
    # cosine distance between centers, based on AS of reads against each center
    vec = get_feature_vector(df_align[['query_id', 'database_id', 'AS']])
    cols = vec.columns[vec.columns != 'id']
    x = vec[cols].values
    cosim = pd.DataFrame(dist_cosine(x.T, x.T), columns=cols)
    cosim['id'] = cols
    df_clst = cluster_optics(cosim, min_samples=2, max_eps=th_m)
    c = (df_clst['cluster_id'] == -1).values
    df_c = df_c[df_c['id'].isin(df_clst['id'].values[c])]
    mqueue = [list(df_clst[df_clst['cluster_id'] == cid]['id'].values)
              for cid in np.unique(df_clst['cluster_id'].values[~c])]
    df_align['weight'] = cluster_weight(df_align, metric='AS')
    logger.info('cluster_merge: %d/%d clusters to merge', np.sum(~c), len(df_clst))
    for i, mout in enumerate(mqueue):
        logger.info('cluster_merge: merging %d clusters, %d/%d', len(mout), i, len(mqueue))
        qlist = cluster_subsample(df_align, mout, N=N, mode='sorted', metric='AS')
        qlist = df_q[df_q['id'].isin(qlist)]
        cout = cluster_compute(qlist, csize, pw_config=pw_config,
                               msa_config=msa_config, tools=tools)
        cout['split'] = True
        df_c = _concat([df_c, cout[['id', 'sequence', 'split']]])
    return _rename_clusters(df_c)


def cluster_weight(df_align, metric='s1'):
    '''Weight of each alignment relative to the read's best alignment.'''
    df_a = get_best(df_align, ['query_id'], metric=metric).rename(columns={metric: 'b1'})
    merged = df_align[['query_id', metric]].merge(df_a[['query_id', 'b1']],
                                                  on='query_id', how='left')
    w = merged[metric].values / (merged[metric].values + merged['b1'].values)
    return w


def cluster_subsample(df_align, dlist, N=100, mode='sorted', metric='AS'):
    '''Pick reads belonging to the given centers, favouring strong alignments.'''
    if len(df_align) == 0:
        return np.array([])
    df_align = df_align.sort_values(by=['database_id', metric])[::-1]  # best first
    x = []
    for d in dlist:
        sub = df_align[df_align['database_id'] == d]
        if mode == 'rand':
            v = sub['query_id'].values
            qlist = v[np.random.permutation(len(v))[:N]]
        else:
            q = sub['query_id'].values
            w = sub['weight'].values.astype(float).copy()
            w[w >= 0.5] = 1
            keep = w > np.random.rand(len(w))  # weighted random acceptance
            q = q[keep]
            if mode == 'semi':
                qlist = q[np.random.permutation(len(q))[:N]]
            else:  # sorted
                qlist = q[:N]
        x += list(qlist)
    return np.unique(x)


def cluster_compute(df_q, csize=20, metric='similarity', pw_config='-k15 -w10 -D',
                    msa_config='-l 0 -r 0', tools=None):
    '''
    Cluster a set of reads: all-vs-all minimap2 -> distance matrix -> OPTICS ->
    spoa consensus per cluster. Returns [id, sequence] of the new centers.
    '''
    empty = pd.DataFrame([], columns=['id', 'sequence'])
    df_q = df_q.copy()
    df_q['id'] = df_q['id'].astype(str)
    if len(df_q) < 2:
        logger.warning('cluster_compute: sample size=%d, too small to do anything', len(df_q))
        return empty
    logger.info('cluster_compute: computing pairwise distance matrix (%d reads)', len(df_q))
    df = run_minimap2(df_q, df_q, config=pw_config + ' --dual=no --for-only',
                      cigar=True, tools=tools)
    if len(df) == 0:
        logger.warning('cluster_compute: not enough samples to cluster')
        return empty
    df = get_best(df, ['query_id', 'database_id'], 'AS')
    df = get_feature_vector(df[['query_id', 'database_id', metric]], symmetric=True)
    df = get_symmetric_matrix(df)
    # similarity -> distance
    cols = df.columns[df.columns != 'id']
    dist = 1 - df[cols].values.astype(float)
    np.fill_diagonal(dist, 0)
    dist_df = pd.DataFrame(dist, columns=cols)
    dist_df['id'] = df['id'].values

    logger.info('cluster_compute: running optics')
    df_clst = cluster_optics(dist_df)
    df_clst = df_clst.merge(df_q, on='id', how='left')
    df_i = df_clst[df_clst['cluster_id'] > -1]
    if len(df_i) == 0:
        logger.warning('cluster_compute: no clusters found')
        return empty
    din = [df_i[df_i['cluster_id'] == cid].iloc[:csize] for cid in np.unique(df_i['cluster_id'])]
    return cluster_spoa_merge(pd.concat(din), config=msa_config, tools=tools)


def cluster_spoa_merge(df, config='-l 0 -r 0', tools=None):
    '''One spoa consensus per cluster_id. Returns [id, sequence].'''
    work = tempfile.mkdtemp(prefix='spoa_', dir=tools['workspace'])
    try:
        df = df.sort_values(by=['cluster_id', 'id'])
        infiles = []
        for cid, grp in df.groupby('cluster_id', sort=True):
            fname = os.path.join(work, str(cid) + '.fa')
            write_fasta(fname, grp[['id', 'sequence']].values)
            infiles.append((cid, fname))
        # run spoa jobs, `threads` at a time
        cmd_base = [tools['spoa']] + config.split()
        nproc = max(1, tools['threads'])
        for start in range(0, len(infiles), nproc):
            procs = []
            for cid, fname in infiles[start:start + nproc]:
                fout = open(fname + '.out', 'w')
                procs.append((subprocess.Popen(cmd_base + [fname], stdout=fout,
                                               stderr=subprocess.PIPE), fout, fname))
            for p, fout, fname in procs:
                _, err = p.communicate()
                fout.close()
                if p.returncode != 0:
                    raise RuntimeError('spoa failed on %s: %s' % (fname, err.decode(errors='replace')))
        data = []
        for cid, fname in infiles:
            seq = read_spoa(fname + '.out')
            if len(seq) == 0:
                # the original kept empty consensus sequences, which later show
                # up in BLAST as "Sequence contains no data"; drop them instead
                logger.warning('cluster_spoa_merge: empty consensus for cluster %s, dropped', cid)
                continue
            data.append(['cluster' + str(cid), seq])
        return pd.DataFrame(data, columns=['id', 'sequence'])
    finally:
        shutil.rmtree(work, ignore_errors=True)


#############################################################################
# Clustering primitives (port of bilge_pype.py)
#############################################################################

def get_best(df, col, metric='AS', stat='idxmax'):
    '''Row with the best `metric` for each group in `col`.'''
    df = df.reset_index(drop=True)
    idx = df.groupby(by=col)[metric].agg(stat).values
    return df.iloc[idx].reset_index(drop=True)


def get_feature_vector(df, symmetric=False):
    '''
    Pivot [query_id, database_id, value] into a matrix.
    returns dataframe with one column per database_id (sorted) plus 'id'
    '''
    q = df.iloc[:, 0].values
    d = df.iloc[:, 1].values
    m = df.iloc[:, 2].values.astype(float)
    if symmetric:
        qlist = np.unique(np.concatenate([q, d]))
        dlist = qlist
    else:
        qlist = np.unique(q)
        dlist = np.unique(d)
    q_dct = {k: i for i, k in enumerate(qlist)}
    d_dct = {k: i for i, k in enumerate(dlist)}
    v = np.zeros([len(qlist), len(dlist)])
    for qi, di, mi in zip(q, d, m):
        v[q_dct[qi], d_dct[di]] = mi
    out = pd.DataFrame(v, columns=list(dlist))
    out['id'] = list(qlist)
    return out


def get_symmetric_matrix(df):
    '''
    Make a square similarity matrix symmetric. For each pair the smaller
    non-zero value is kept (0 = no alignment), as in bilge_pype with
    sym_larger=False.
    '''
    cols = df.columns[df.columns != 'id']
    a = df[cols].values.astype(float)
    at = a.T
    sym = np.where(a == 0, at, np.where(at == 0, a, np.minimum(a, at)))
    out = pd.DataFrame(sym, columns=cols)
    out['id'] = df['id'].values
    return out


def dist_cosine(v1, v2):
    '''Cosine distance between rows of v1 and rows of v2 (v1 == v2 here).'''
    B = np.matmul(v1, v2.T)
    norms = np.sqrt(np.diag(B)).reshape(-1, 1)
    D = B / np.matmul(norms, norms.T)
    D[D > 1] = 1  # numerical imprecision
    return 1 - D


def cluster_optics(df, min_samples=None, max_eps=None, xi=0.05,
                   cluster_method='dbscan', min_cluster_size=None):
    '''
    OPTICS on a precomputed distance matrix (df = [d_1..d_n, id]).
    If min_samples is None, search for the value giving the most clusters with
    the fewest outliers (same search as bilge_pype.cluster_OPTICS).
    returns dataframe [id, cluster_id, ordering, reachability] in OPTICS order
    '''
    if len(df) < 3:
        logger.warning('cluster_optics: sample size too small, each item is its own cluster')
        return pd.DataFrame({'id': df['id'].values, 'cluster_id': np.arange(len(df)),
                             'ordering': np.arange(len(df)), 'reachability': -1.0})
    cols = df.columns[df.columns != 'id']
    X = df[cols].values.astype(float)
    if max_eps is None:
        z = np.abs(X)
        max_eps = (np.max(z) + np.min(z)) / 2
    if min_samples is None or min_samples > len(df):
        min_samples = np.ceil(len(df) / 2)
        max_iter = 100
    else:
        max_iter = 1

    prev_nout, prev_nclust, prev_min = len(df), 0, min_samples * 2
    best_clust, best_nclust, best_nout = None, 0, prev_nout
    delta = min_samples / 2
    tmp = prev_min
    for i in range(max_iter):
        min_samples = int(min(min_samples, len(df)))  # sklearn needs min_samples <= n
        clust = sklearn.cluster.OPTICS(min_samples=min_samples, min_cluster_size=min_cluster_size,
                                       xi=xi, max_eps=max_eps, cluster_method=cluster_method,
                                       metric='precomputed')
        with warnings.catch_warnings():
            # the min_samples search deliberately tries settings where every
            # point is an outlier; sklearn warns about those, which is expected
            warnings.simplefilter('ignore', UserWarning)
            clust.fit(X)
        nout = int(np.sum(clust.labels_ == -1))
        nclust = int(np.max(clust.labels_) + 1)
        logger.debug('cluster_optics: iter=%d min_samples=%d clusters=%d outliers=%d',
                     i, min_samples, nclust, nout)
        if nclust > prev_nclust or (nclust == prev_nclust and nout <= prev_nout):
            tmp = min_samples
            delta = int((prev_min - min_samples) / 2)
            min_samples -= delta
        elif nclust < prev_nclust:
            tmp = min_samples
            delta = int((prev_min + min_samples) / 2) - min_samples
            min_samples = prev_min + delta
        if best_clust is None or nclust > best_nclust or (nclust == best_nclust and nout <= best_nout):
            if np.min(clust.reachability_) >= 0:
                best_clust, best_nclust, best_nout = clust, nclust, nout
        if delta == 0 or min_samples <= 2 or nout == 0:
            break
        prev_nclust, prev_nout, prev_min = nclust, nout, tmp

    s = best_clust.ordering_
    reach = best_clust.reachability_[s].copy()
    finite = reach[np.isfinite(reach)]
    reach[~np.isfinite(reach)] = np.max(finite) if len(finite) else -1
    return pd.DataFrame({'id': df['id'].values[s], 'cluster_id': best_clust.labels_[s].astype(int),
                         'ordering': np.arange(len(s)), 'reachability': reach.astype(float)})


def cluster_kmeans_1d(values, n_clusters, seed=None):
    '''
    k-means on a 1-D array; labels are re-numbered so 0 is the partition with
    the lowest mean value (bilge_pype.cluster_Kmeans with ordered=True).
    '''
    x = values.reshape(-1, 1)
    km = sklearn.cluster.MiniBatchKMeans(n_clusters=n_clusters, n_init=10, max_iter=100,
                                         random_state=seed)
    labels = km.fit_predict(x)
    present = np.unique(labels)
    means = [np.mean(values[labels == l]) for l in present]
    order = present[np.argsort(means)]
    new = np.zeros(len(labels), dtype=int)
    for rank, l in enumerate(order):
        new[labels == l] = rank
    return new


#############################################################################
# minimap2 / spoa / file helpers
#############################################################################

PAF_COLS = ['query_id', 'q_len', 'q_start', 'q_end', 'orientation', 'database_id',
            't_len', 't_start', 't_end', 'match', 'tot', 'mapq']
PAF_TAGS = ['tp', 'cm', 's1', 's2', 'NM', 'AS', 'ms', 'nn', 'rl', 'cg']
PAF_INT = ['q_len', 'q_start', 'q_end', 't_len', 't_start', 't_end', 'match', 'tot',
           'mapq', 'cm', 's1', 's2', 'NM', 'AS', 'ms', 'nn', 'rl']


def run_minimap2(query, database, config='-x map-ont', cigar=True, tools=None):
    '''
    Align query [id, sequence] against database [id, sequence] with minimap2.
    returns dataframe of PAF records plus match_score and similarity
    '''
    out_cols = PAF_COLS + PAF_TAGS[:-1] + ['CIGAR', 'match_score', 'similarity']
    if len(query) == 0 or len(database) == 0:
        return pd.DataFrame([], columns=out_cols)
    work = tempfile.mkdtemp(prefix='mm2_', dir=tools['workspace'])
    try:
        # minimap2 truncates ids at whitespace, so use safe temporary ids.
        # The same id must get the same name in query and database, because
        # the -D / --dual=no options compare query and target names.
        qmap = _safe_names(query['id'].values)
        dmap = _safe_names(database['id'].values)
        qfile = os.path.join(work, 'query.fa')
        dbfile = os.path.join(work, 'database.fa')
        outfile = os.path.join(work, 'results.paf')
        write_fasta(qfile, np.column_stack([list(qmap.keys()), query['sequence'].values]))
        write_fasta(dbfile, np.column_stack([list(dmap.keys()), database['sequence'].values]))
        cmd = [tools['minimap2'], '-t', str(tools['threads'])]
        if cigar:
            cmd.append('-c')
        cmd += config.split() + [dbfile, qfile, '-o', outfile]
        res = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        if res.returncode != 0:
            raise RuntimeError('minimap2 failed: ' + res.stderr.decode(errors='replace')[-2000:])
        data = parse_paf(outfile)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    if len(data) == 0:
        return pd.DataFrame([], columns=out_cols)
    data['query_id'] = data['query_id'].map(qmap)
    data['database_id'] = data['database_id'].map(dmap)
    for c in PAF_INT:
        data[c] = data[c].astype(int)
    data['match_score'] = data['match'] / data['tot']
    shorter = np.minimum(data['q_len'].values, data['t_len'].values)
    data['similarity'] = np.clip(1 - data['NM'].values / shorter, 0, None)
    if not cigar:
        data = data.drop(columns=['CIGAR'])
    return data


def _safe_names(ids):
    '''map whitespace-free name -> original id (index-prefixed if not unique)'''
    names = [re.sub(r'\s+', '_', str(i)) for i in ids]
    if len(set(names)) != len(names):
        names = [str(k) + '_' + n for k, n in enumerate(names)]
    return dict(zip(names, ids))


def parse_paf(fname):
    '''Parse a minimap2 PAF file into a dataframe (missing tags -> 0).'''
    key = {t: i for i, t in enumerate(PAF_TAGS)}
    rows = []
    with open(fname) as f:
        for line in f:
            out = line.rstrip('\n').split('\t')
            if len(out) < len(PAF_COLS):
                continue
            info = [0] * len(PAF_TAGS)
            for field in out[len(PAF_COLS):]:
                parts = field.split(':')
                if parts[0] in key:
                    info[key[parts[0]]] = parts[-1]
            rows.append(out[:len(PAF_COLS)] + info)
    df = pd.DataFrame(rows, columns=PAF_COLS + PAF_TAGS)
    return df.rename(columns={'cg': 'CIGAR'})


def read_spoa(fname):
    '''Consensus = last record of spoa output, gaps removed.'''
    recs = read_fasta(fname)
    if len(recs) == 0:
        return ''
    return recs['sequence'].values[-1].replace('-', '')


def write_fasta(fname, data):
    '''data = array of [id, sequence] rows'''
    with open(fname, 'w') as f:
        for rid, seq in data:
            f.write('>' + str(rid) + '\n' + str(seq) + '\n')


def read_fasta(fname):
    '''Read a fasta file into [id, sequence].'''
    ids, seqs, cur = [], [], []
    opener = gzip.open if fname.endswith('.gz') else open
    with opener(fname, 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if ids:
                    seqs.append(''.join(cur))
                ids.append(line[1:])
                cur = []
            elif line:
                cur.append(line.strip())
    if ids:
        seqs.append(''.join(cur))
    return pd.DataFrame({'id': ids, 'sequence': seqs})


def read_fastq(fname):
    '''Read a fastq file into [id, sequence] (qualities are not needed here).'''
    ids, seqs = [], []
    opener = gzip.open if fname.endswith('.gz') else open
    with opener(fname, 'rt') as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            if header.startswith('@'):
                ids.append(header[1:].strip())
                seqs.append(seq)
    return pd.DataFrame({'id': ids, 'sequence': seqs})


def _load_sequences(fname):
    if not os.path.isfile(fname):
        raise FileNotFoundError(fname)
    low = fname.lower()
    if re.search(r'\.(fa|fasta|fna)(\.gz)?$', low):
        df = read_fasta(fname)
    elif re.search(r'\.(fq|fastq)(\.gz)?$', low):
        df = read_fastq(fname)
    else:
        df = pd.read_csv(fname)
    if 'sequence' not in df.columns:
        raise ValueError(fname + ' needs a "sequence" column')
    df = df[df['sequence'].notna()].copy()
    df['sequence'] = df['sequence'].astype(str).str.strip().str.upper()
    df = df[df['sequence'].str.len() > 0]
    if 'id' not in df.columns:
        df['id'] = ['seq' + str(i) for i in range(len(df))]
    return df.reset_index(drop=True)


def _concat(frames):
    '''pd.concat that skips empty frames (avoids pandas FutureWarnings)'''
    keep = [f for f in frames if len(f) > 0]
    if not keep:
        return frames[0].iloc[0:0]
    return pd.concat(keep)


def _rename_clusters(df_c):
    df_c = df_c.reset_index(drop=True)
    df_c['id'] = ['cluster' + str(i) for i in range(len(df_c))]
    return df_c


def _check_tool(tool):
    if shutil.which(tool) is None:
        raise FileNotFoundError('%s not found on PATH; install it or pass its full path' % tool)


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description='Cluster amplicon reads into consensus sequences')
    p.add_argument('-i', dest='input_file', required=True)
    p.add_argument('-o', dest='output_file', required=True)
    p.add_argument('-iter', dest='n_iter', type=int, default=DEFAULTS['n_iter'])
    p.add_argument('-seed', dest='seed', type=int, default=None)
    p.add_argument('-t', dest='threads', type=int, default=3)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')
    cluster(a.input_file, a.output_file, n_iter=a.n_iter, seed=a.seed, threads=a.threads)
