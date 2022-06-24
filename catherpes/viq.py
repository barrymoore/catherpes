#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created on Tue Apr 26 21:05:28 UTC 2022

Synopsis:

viq.py viq_output.txt

Description:

This python module is intended to be imported, however it can be
tested by passing a vIQ output file, in which case it will print the
#A section formated as TSV to STDOUT.

Positional Arguments:

  file: A vIQ output file.

"""

import sys
import re
import json
import collections
import pandas as pd

class VIQ():

    gene_keys = ['rank', 'chr', 'gene', 'transcript', 'vid',
                 'csq', 'dist', 'denovo', 'type', 'zygo', 'csn',
                 'pldy', 'sites', 'par', 'loc', 'length', 'gqs',
                 'gflg', 'gflpr', 'ppp', 'vpene', 'breath', 'fix',
                 'viqscr', 'p_scor', 's_scor', 'phev_k',
                 'vvp_svp', 'vaast', 'rprob', 'g_tag', 'p_mod',
                 's_mod', 'g_tag_scr', 'clinvar', 'var_qual',
                 'rid', 'loc2', 'maf', 'incndtl', 'payload']

    mim_keys = ['rank', 'mim', 'gene', 'phev_gene', 'phev_mim',
                'mim_mpr', 'imprm', 'viqscr', 'miqscr', 'inc', 'disease']

    multi_keys = ['rank', 'mim', 'vid', 'num_g', 'type', 'zyg',
                  'csqs', 'ploidy', 'phev_gene', 'phev_mim', 'mim_mpr',
                  'viqscr', 'miqscr', 'inc', 'chr', 'name', 'best_mim', 'genes']

    def __init__(self, file=None):

        # Function to create dict autovivification - come on Python, really?
        def auto_dict():
            return collections.defaultdict(auto_dict)

        self.file = file
        self.gene_records = dict()
        self.mim_records = dict()
        self.multi_records = dict()
        self.hpo = list()
        self.meta = {}
        self.meta['sv_prior_obs'] = auto_dict()
        self.meta['ave_depth_of_coverage'] = auto_dict()
        self.meta['eof'] = False
        self.parse_file()

    def parse_file(self):

        # Rank CHR Gene Transcript vID CSQ DIST Denovo Type Zygo CSN PLDY
        # SITES Par Loc Length GQS GFLG GFLpr PPP vPene breath FIX vIQscr
        # p_scor s_scor PHEV/K VVP/SVP VAAST RPROB G_tag p_mod s_mod
        # G_tag_scr ClinVar var_qual vID CHR;BEG;END

        GeneRecord =collections.namedtuple('GeneRecord',
                                ' '.join(VIQ.gene_keys))

        MimRecord =collections.namedtuple('MimRecord',
                                ' '.join(VIQ.mim_keys))

        MultiRecord =collections.namedtuple('MultiRecord',
                                ' '.join(VIQ.multi_keys))

        with open(self.file, "r") as fh:
            for line in fh:
                line = line.strip()

                if re.search('^\#\#', line):
                    pass
                if re.search('#A\tRank\tCHR\tGene', line):
                    pass
                elif re.match("#A", line):
                    a_values = [x.strip() for x in line.split('\t')]
                    # Get rid of the leading #A column
                    a_values.pop(0)
                    col_count = len(a_values)
                    if (col_count != 38 & col_count != 39):
                        sys.exit('FATAL : incorrect_column_count_section_A, (expected 38 ' +
                                 f'or 39, but got {col_count} columns) {line}')

                    gene = a_values[2]

                    rec_dict = dict(zip(VIQ.gene_keys, a_values[0:38]))

                    if col_count == 39:
                        rec_dict['payload'] = a_values[38]
                    else:
                        rec_dict['payload'] = None

                    # Parse denovo
                    (rec_dict['denovo'], rec_dict['maf']) = rec_dict['denovo'].split('(')
                    rec_dict['maf'] = re.sub('\)$', '', rec_dict['maf'])

                    # Parse indendental
                    rec_dict['incndtl'] = None
                    imatch = re.search('[gpn]$', rec_dict['clinvar'])
                    if imatch:
                        rec_dict['clinvar'] = re.sub('([gpn])$', '', rec_dict['clinvar'])
                        rec_dict['incndtl'] = imatch.group(1)

                    # Parse var_qual
                    rec_dict['var_qual'] = re.sub('^\(', '', rec_dict['var_qual'])
                    vq_data = []
                    vq_dict = dict()
                    if re.search('\s', rec_dict['var_qual']):
                        vq_data = rec_dict['var_qual'].split()
                        vq_dict['type'] = 'sv'
                        vq_dict['values'] = vq_data.pop(0)
                        for vq_pair in vq_data:
                            (key, value) = vq_pair.split(':')
                            vq_dict[key] = value
                    else:
                        vq_data = rec_dict['var_qual'].split('|')
                        vq_data[0] = vq_data[0].split(':')
                        vq_dict['type'] = 'snv'
                        vq_dict['values'] = vq_data

                    rec_dict['var_qual'] = vq_dict

                    # Parse g_tag_scr: '0.5;0.5'
                    rec_dict['g_tag_scr'] = rec_dict['g_tag_scr'].split(';')

                    # Parse 'loc2': '20;10616332;10656694'
                    rec_dict['loc2'] = dict(zip(['chr','start','end'], rec_dict['loc2'].split(';')))

                    # Parse 'pldy': '1(0.90195 0.92688)'
                    (my_pldy, the_rest) = rec_dict['pldy'].split('(')
                    the_rest = re.sub('\)', '', the_rest)
                    pldy_values = the_rest.split()
                    rec_dict['pldy'] = dict([('pldy', my_pldy), ('values', pldy_values)]) 
                    
                    record = GeneRecord(**rec_dict)
                    self.gene_records[gene] = record

                elif re.match("#B", line):
                    ##------------------------------------------------------------------------------
                    ##------------------------ POSSIBLE MENDELIAN DIAGNOSES ------------------------
                    ##------------------------------------------------------------------------------
                    ## RANK  MIM     GENE     PHEV_GENE PHEV_MIM  MIM_MPR IMPRM   vIQscr  mIQscr  INC  DISEASE
                    #B 0     118450  JAG1     0.957     0.992     0.987   0.5     1.126   1.545   N    AWS
                    #B 1     618786  SUZ12    0.952     0.904     0.9     0.5     1.21    0.945   N    Imagawa-Matsumoto syndrome
                    #B 2     617140  SON      0.953     0.979     0.975   0.5     -0.445  0.212   N    ZTTK SYNDROME
                    #B 3     617159  CHD4     0.878     0.836     0.832   0.5     0.225   0.055   N    SIFRIM-HITZ-WEISS SYNDROME

                    b_values = [x.strip() for x in line.split('\t')]
                    # Get rid of the leading #B column
                    b_values.pop(0)
                    col_count = len(b_values)
                    if col_count != 11:
                        sys.exit('FATAL : incorrect_column_count_section_B, (expected 11 ' +
                                 f', but got {col_count} columns) {line}')

                    gene = b_values[2]

                    rec_dict = dict(zip(VIQ.mim_keys, b_values))

                    record = MimRecord(**rec_dict)
                    self.mim_records[gene] = record

                elif re.match("#C", line):
                    ##------------------------------------------------------------------------------
                    ##------------------------ POSSIBLE MULTIGENIC DIAGNOSES -----------------------
                    ##------------------------------------------------------------------------------
                    ##  RANK  MIM     VID              NUM_G  TYPE  ZYG  CSQS  PLOIDY  PHEV_GENE  PHEV_MIM  MIM_MPR vIQscr  mIQscr  INC  CHR
                    #C  0     113100  BadgeGrp37362.2  2      4     1    1     1       0.925      0.662     0.66    -1.291  -1.319  N    20
                    #C  1     NO_MIM  BadgeGrp7940.2   2      4     1    1     1       0.247      0.012     0.002   -1.612  -3.12   Y    16
                    #C  2     NO_MIM  BadgeGrp6437.2   2      4     1    1     1       0.297      0.001     0.002   -1.722  -3.199  Y    2
                    #C  3     NO_MIM  BadgeGrp2507.3   3      6     1    8     3       0.264      0         0.002   -1.968  -3.375  Y    7
                    c_values = [x.strip() for x in line.split('\t')]
                    # Get rid of the leading #B column
                    c_values.pop(0)
                    col_count = len(c_values)
                    if col_count != 18:
                        sys.exit('FATAL : incorrect_column_count_section_C, (expected 18 ' +
                                 f', but got {col_count} columns) {line}')

                    vid = c_values[2]

                    rec_dict = dict(zip(VIQ.multi_keys, c_values))

                    record = MultiRecord(**rec_dict)
                    self.multi_records[vid] = record

                ##  0  HP:0011121  1
                elif re.search("^##\s+\d+\s+HP:\d+", line):
                    m = re.search("^##\s+\d+\s+(HP:\d+)", line)
                    self.hpo.append(m.group(1))

                # elif re.search("^XXX", line):
                # pass
                # 	     ## GENOTYPE SKEW CHECK. P_value alpha = 0.00217391 based upon 90 previous observations.
                # 	     ## CHR  	NUM_HET 	NUM_HOM 	%HOM_OBS	%HOM_EXP	P_VALUE      	 het <- SKEW -> hom
                # 	     ## 1    	989     	475     	0.307122	0.274696	0.107468	          |
                # 	     ## 10   	279     	172     	0.307092	0.291090	0.317974	          |
                # 	     ## 11   	460     	360     	0.361816	0.301593	0.103433	          |
                # 	     ## 12   	397     	256     	0.258878	0.303959	0.089466	          |
                # 	     ## 13   	154     	83      	0.223900	0.263517	0.209185	          |
                # 	     ## 14   	221     	136     	0.257753	0.303610	0.190588	          |
                # 	     ## 15   	269     	148     	0.242351	0.311576	0.071857	          |
                # 	     ## 16   	245     	174     	0.383414	0.260717	0.000274	          |+++
                # 	     ## 17   	451     	273     	0.268743	0.297246	0.197559	          |
                # 	     ## 18   	101     	96      	0.381859	0.293539	0.043265	          |
                # 	     ## 19   	672     	386     	0.287421	0.292273	0.439947	          |
                # 	     ## 2    	511     	320     	0.327274	0.293646	0.116324	          |
                # 	     ## 20   	212     	89      	0.233125	0.289368	0.104478	          |
                # 	     ## 21   	134     	43      	0.137459	0.288819	0.003829	          |
                # 	     ## 22   	171     	93      	0.260173	0.306226	0.184831	          |
                # 	     ## 3    	407     	255     	0.283468	0.297334	0.353994	          |
                # 	     ## 4    	277     	232     	0.444243	0.307158	0.000250	          |+++
                # 	     ## 5    	310     	225     	0.349154	0.303541	0.128416	          |
                # 	     ## 6    	644     	323     	0.254281	0.281636	0.317448	          |
                # 	     ## 7    	433     	173     	0.199651	0.280076	0.019431	          |
                # 	     ## 8    	227     	152     	0.290023	0.295137	0.446953	          |
                # 	     ## 9    	263     	184     	0.295424	0.292450	0.468794	          |
                # 	     ## X    	20      	126     	0.872348	0.647167	0.223214	          |
                # 	     ##

                ## LOH DETECTED:YES
                elif re.match("## LOH DETECTED", line):
                    m = re.search("^## LOH DETECTED:(.*)", line)
                    self.meta['loh_detected'] = m.group(1)

                ## SKEW DETECTED:YES (7.16236527310748)
                elif re.match("^## SKEW DETECTED:", line):
                    m = re.search("^## SKEW DETECTED:(\S+)\s+(\(.*\))", line)
                    self.meta['skew_detected'] = m.group(1)
                    self.meta['skew'] = m.group(2)

                ## Parental relationships MSUM:1 FSUM:1 FMR:1 DSUM:778 DSR:389
                elif re.match("## Parental relationships MSUM", line):
                    m = re.search("^## Parental relationships (MSUM:(\d+)\s+FSUM:(\d+)\s+FMR:(\S+)\s+DSUM:(\d+)\s+DSR:(\d+))", line)
                    self.meta['parental_relationships'] = m.group(1)
                    self.meta['pr_msum'] = m.group(2)
                    self.meta['pr_fsum'] = m.group(3)
                    self.meta['pr_fmr']  = m.group(4)
                    self.meta['pr_dsum'] = m.group(5)
                    self.meta['pr_dsr']  = m.group(6)

                ## Estimated Consanguinity RAW:5.57 % Ancestry ADJ:0 %
                elif re.match("## Estimated Consanguinity", line):
                    m = re.search("^## Estimated Consanguinity\s+(RAW:(.*?)\s+%\s+Ancestry\s+ADJ:(.*?)\s+%)", line)
                    self.meta['estimated_consanguinity'] = m.group(1)
                    self.meta['ec_raw'] = m.group(2)
                    self.meta['ec_adj'] = m.group(3)

                ## CSN PRIOR:0.5
                elif re.match("## CSN PRIOR", line):
                    m = re.search("^## CSN PRIOR:(.*)", line)
                    self.meta['csn_prior'] = m.group(1)

                ## VARIANTS_IN:14993 NUMBER OF SVs (rows):1754 PASSING_FILTERS:2307 p_obs:0.5
                elif re.match("## VARIANTS_IN", line):
                    m = re.search("^## VARIANTS_IN:((\d+)\s+NUMBER\s+OF\s+SVs\s+\(rows\):(\d+)\s+PASSING_FILTERS:(\d+)\s+p_obs:(.+))", line)
                    self.meta['variants_in'] = m.group(1)
                    self.meta['vars_snv'] = m.group(2)
                    self.meta['vars_sv'] = m.group(3)
                    self.meta['vars_passing'] = m.group(4)
                    self.meta['vars_p_obs'] = m.group(5)

                ## Ext. SVs in file:1657
                elif re.match("## Ext. SVs in file", line):
                    m = re.search("^## Ext. SVs in file:(\d+)", line)
                    self.meta['ext_svs'] = m.group(1)

                ## Ext. SVs in kept:0
                elif re.match("## Ext. SVs in kept", line):
                    m = re.search("^## Ext. SVs in kept:(\d+)", line)
                    self.meta['ext_svs_kept'] = m.group(1)

                ## Badges trimmed:0
                elif re.match("## Badges trimmed", line):
                    m = re.search("^## Badges trimmed:(\d+)", line)
                    self.meta['badges_trimmed'] = m.group(1)

                ## Badges considered:29788
                elif re.match("## Badges considered", line):
                    m = re.search("^## Badges considered:(\d+)", line)
                    self.meta['badges_considered'] = m.group(1)

                ## Badges kept:5
                elif re.match("## Badges kept", line):
                    m = re.search("^## Badges kept:(\d+)", line)
                    self.meta['badges_kept'] = m.group(1)

                ## Badges deleted:1
                elif re.match("## Badges deleted", line):
                    m = re.search("^## Badges deleted:(\d+)", line)
                    self.meta['badges_deleted'] = m.group(1)

                ## Proband Ancestry:Other
                ## Relative -Log Likelihoods:
                ## Other:10355,European(non-Finish):11019,Finish:11021,Ashkenazi:11079,Asian:11373,African:12613
                elif re.match("## Proband Ancestry", line):
                    m = re.search("^## Proband Ancestry:(.*?)\s+Relative\s+-Log\s+Likelihoods:\s+(.*)", line)
                    self.meta['proband_ancestry'] = m.group(1)
                    pa_rrl_txt = m.group(2)
                    pa_rll_list = pa_rrl_txt.split(',')
                    pa_rll_dict = dict()
                    for pa_rll in pa_rll_list:
                        (ans, rll) = pa_rll.split(':')
                        pa_rll_dict[ans] = rll
                    self.meta['pa_rlls'] = pa_rll_dict

                ## Proband Sex:f P(MALE):0.322425452976684
                elif re.match("## Proband Sex", line):
                    m = re.search("^## Proband Sex:([fm])\s+P\(MALE\):(.*)", line)
                    self.meta['proband_sex'] = m.group(1)
                    self.meta['ps_prob_male'] = m.group(2)

                ## NUM HPO TERMS:3
                elif re.match("## NUM HPO TERMS", line):
                    m = re.search("^## NUM HPO TERMS:(\d+)", line)
                    self.meta['hpo_term_count'] = m.group(1)

                ## ADJ FOR INBREEDING:0.541174032480132
                elif re.match("## ADJ FOR INBREEDING", line):
                    m = re.search("^## ADJ FOR INBREEDING:(.*)", line)
                    self.meta['adj_for_inbreeding'] = m.group(1)

                ## MODE:SINGLETON
                elif re.match("## MODE", line):
                    m = re.search("^## MODE:(.*)", line)
                    self.meta['mode'] = m.group(1)

                ## Number of variants failing -e m  a:cov:120,bias:171,tot:248 x:cov:5,bias:9,tot:14
                elif re.match("## Number of variants failing -e m", line):
                    m = re.search("^## Number of variants failing -e m  (.*)", line)
                    self.meta['num_vars_failing_em'] = m.group(1)
                    data = self.meta['num_vars_failing_em'].split()
                    this_type = ''
                    for item in data:
                        m = re.match('^([ax]):', item)
                        this_type = m.group(1)
                        
                        item = re.sub('^[ax]:', '', item)

                        nvf_dict = dict()
                        for e in item.split(','):
                            (key, value) = e.split(':')
                            nvf_dict[key] = value
                            self.meta['nvf_' + this_type] = nvf_dict

                ## VAAST-VVP COOR:0.51
                elif re.match("## VAAST-VVP COOR", line):
                    m = re.search("^## VAAST-VVP COOR:(.*)", line)
                    self.meta['vaast_vvp_coor'] = m.group(1)

                ## BLS-BND COOR:0.0049079754601227
                elif re.match("## BLS-BND COOR", line):
                    m = re.search("^## BLS-BND COOR:(.*)", line)
                    self.meta['bls_bnd_coor'] = m.group(1)

                ## BLS-NOA COOR:0.0049079754601227
                elif re.match("## BLS-NOA COOR", line):
                    m = re.search("^## BLS-NOA COOR:(.*)", line)
                    self.meta['bls_noa_coor'] = m.group(1)

                ## PHEV-KPR COOR:0.34
                elif re.match("## PHEV-KPR COOR", line):
                    m = re.search("^## PHEV-KPR COOR:(.*)", line)
                    self.meta['phev_kpr_coor'] = m.group(1)

                ## CLIN-VVP-VAAST:0.725758988640042
                elif re.match("## CLIN-VVP-VAAST", line):
                    m = re.search("^## CLIN-VVP-VAAST:(.*)", line)
                    self.meta['clin_vvp_vaast'] = m.group(1)

                ## COVERAGE-HOMOZYGOSITY JACC:0
                elif re.match("## COVERAGE-HOMOZYGOSITY JACC", line):
                    m = re.search("^## COVERAGE-HOMOZYGOSITY JACC:(.*)", line)
                    self.meta['coverage_homozygosity_jacc'] = m.group(1)

                ## COVERAGE-HETEROZYGOSITY JACC:0.04
                elif re.match("## COVERAGE-HETEROZYGOSITY JACC", line):
                    m = re.search("^## COVERAGE-HETEROZYGOSITY JACC:(.*)", line)
                    self.meta['coverage_heterozygosity_jacc'] = m.group(1)

                ## PHEV-GENE-MIM COOR:0.42
                elif re.match("## PHEV-GENE-MIM COOR", line):
                    m = re.search("^## PHEV-GENE-MIM COOR:(.*)", line)
                    self.meta['phev_gene_mim_coor'] = m.group(1)

                ## K_PRIOR:0.8344
                elif re.match("## K_PRIOR:", line):
                    m = re.search("^## K_PRIOR:(.*)", line)
                    self.meta['k_prior'] = m.group(1)

                ## MIM-GEN SCORES COOR:0.08
                elif re.match("## MIM-GEN SCORES COOR", line):
                    m = re.search("^## MIM-GEN SCORES COOR:(.*)", line)
                    self.meta['mim_gen_scores_coor'] = m.group(1)

                ## K_PRIOR_PROB:0.232164449818622
                elif re.match("## K_PRIOR_PROB", line):
                    m = re.search("^## K_PRIOR_PROB:(.*)", line)
                    self.meta['k_prior_prob'] = m.group(1)

                ## U_PRIOR_PROB:0.405990586221652
                elif re.match("## U_PRIOR_PROB", line):
                    m = re.search("^## U_PRIOR_PROB:(.*)", line)
                    self.meta['u_prior_prob'] = m.group(1)

                # SV PRIOR OBS	a	f	1	80
                elif re.match("# SV PRIOR OBS", line):
                    m = re.search("^# SV PRIOR OBS\s+([ax])\s+([fm])\s+(\d+)\s+(\d+)", line)

                    if 'sv_prior_obs' not in self.meta:
                        self.meta['sv_prior_obs'] = {}
                    chrom = m.group(1)
                    sex = m.group(2)

                    self.meta['sv_prior_obs'][chrom][sex]['_first'] = m.group(3)
                    self.meta['sv_prior_obs'][chrom][sex]['_second'] = m.group(4)

                ## AVE DEPTH OF COVERAGE 1 MEAN:38.2462006079027 VARIANCE:266.838590703536
                elif re.match("## AVE DEPTH OF COVERAGE", line):
                    m = re.search("^## AVE DEPTH OF COVERAGE\s+(\S+)\s+MEAN:(\S+)\s+VARIANCE:(\S+)", line)
                    chrom = m.group(1)
                    mean = m.group(2)
                    var = m.group(3)

                    self.meta['ave_depth_of_coverage'][chrom]['mean'] = mean
                    self.meta['ave_depth_of_coverage'][chrom]['variance'] = var

                ## CLIN PRIORS
                ##	0	0.00502478219010376
                ##	1	0.00503303631683697
                ##	2	0.00504954457030337
                ##	3	0.00508256107723618
                ##	5	0.994983471936629
                ##	4	0.994983471936629
                ##	9	0.994983471936629
                ##	6	0.994983471936629
                ##	8	0.994983471936629
                ##
                ## CSQS
                ##	12,33	3.35638047929113e-05
                ##	3,23	3.35638047929113e-05
                ##	4,5	3.35638047929113e-05
                ##	5,13	3.35638047929113e-05
                ##	7	3.35638047929113e-05
                ##	2	0.000100691414378734
                ##	13,20	0.000134255219171645
                ##	1	0.000234946633550379
                ##	3	0.000234946633550379
                ##	11,13	0.000268510438343291
                ##	13,17	0.000302074243136202
                ##	4	0.000302074243136202
                ##	10	0.00050345707189367
                ##	9	0.000704839900651138
                ##	12,35	0.000738403705444049
                ##	5	0.00104047794858025
                ##	13,23	0.00288648721219037
                ##	11	0.0192320601463382
                ##

                ## TYPE FREQUENCIES:	1:0.862332695984704	2:0.137667304015296	SUM:523
                elif re.match("## TYPE FREQUENCIES", line):
                    m = re.search("^## TYPE FREQUENCIES:\s+(.*)\s+SUM:(\d+)", line)

                    if 'type_frequencies' not in self.meta:
                        self.meta['type_frequencies'] = {}
                    type_txt = m.group(1)
                    for my_type in type_txt.split():
                        (k, v) = my_type.split(':')
                        self.meta['type_frequencies'][k] = v

                    self.meta['type_frequencies']['sum'] = m.group(2)

                ## TFA (INDEL) ADJUSTMENT	STATUS:on	TFA:0.389454872134458	Pval:0.000211756350740666	exp freq:0.0878152942746563
                elif re.match("## TFA \(INDEL\) ADJUSTMENT", line):
                    m = re.search("^## TFA \(INDEL\) ADJUSTMENT\s+(STATUS:(\S+)\s+TFA:(\S+)\s+Pval:(\S+)\s+exp\s+freq:(\S+))", line)
                    self.meta['tfa_adjustment'] = m.group(1)
                    self.meta['tfa_status'] = m.group(2)
                    self.meta['tfa_tfa'] = m.group(3)
                    self.meta['tfa_pval'] = m.group(4)
                    self.meta['tfa_exp_freq'] = m.group(5)

                ## LOW QUALITY LIST FILE DETECTED:NO. Frac Vars w/out reads:0
                elif re.match("## LOW QUALITY LIST FILE DETECTED", line):
                    m = re.search("^## LOW QUALITY LIST FILE DETECTED:(\S+)\s+Frac Vars w/out reads:(\S+)", line)
                    self.meta['low_quality_list_file'] = m.group(1)
                    self.meta['low_quality_list_file_frac'] = m.group(2)

                ## NULL READ COUNT ADJ. With:1 Without:1
                elif re.match("## NULL READ COUNT ADJ\.", line):
                    m = re.search("^## NULL READ COUNT ADJ\.\s+With:(\S+)\s+Without:(\S+)", line)
                    self.meta['null_read_count_adj_with'] = m.group(1)
                    self.meta['null_read_count_adj_without'] = m.group(2)

                ## CMD:/home/ubuntu/vIQ/bin/vIQ2stripe -a /home/ubuntu/fabric_viq_workflow/snakemake/rady_benchmarks/viq.config -b 0.95 -B viq_ontos/bad.txt -c  -d  -T WGS -e m -f 0.005 -g  -H  -h  -k  -l VIQ/coding_dist.CT.CS4EBE.viq_list.txt -m s -M Phevor/phv_renorm_mim.CT.CS4EBE.txt -o  -p 0.5 -Q m -q n -r n -s Phevor/phv_renorm.CT.CS4EBE.txt -v  -w  -x  -y  -z
                elif re.match("## CMD:", line):
                    m = re.search("^## CMD:(.*)", line)
                    self.meta['cmd'] = m.group(1)

                ## Threshold for SV support (BFT):0.698970004336019
                elif re.match("## Threshold for SV support \(BFT\)", line):
                    m = re.search("^## Threshold for SV support \(BFT\):(\S+)", line)
                    self.meta['threshold_sv_support'] = m.group(1)

                ## K:6.51329346539193
                elif re.match("## K:", line):
                    m = re.search("^## K:(\S+)", line)
                    self.meta['k'] = m.group(1)

                ## BEE:0
                elif re.match("## BEE:", line):
                    m = re.search("^## BEE:(.*)", line)
                    self.meta['bee'] = m.group(1)

                ## VERSION:7.0.02
                elif re.match("## VERSION:", line):
                    m = re.search("^## VERSION:(\S+)", line)
                    self.meta['version'] = m.group(1)

                ## GMT:Tue Jun 21 22:32:36 2022
                elif re.match("## GMT:", line):
                    m = re.search("^## GMT:(.*)", line)
                    self.meta['gmt'] = m.group(1)

                ## EOF
                elif re.match("## EOF", line):
                    m = re.search("^## EOF", line)
                    self.meta['eof'] = True

    def genes_as_list(self):
        df = self.genes_as_df()
        return df.values.tolist()
    
    def genes_as_json(self):
        return json.dumps(self.gene_records, indent=4, sort_keys=True)
    
    def genes_as_df(self):

        genes = []
        for record in self.gene_records.values():
            gene = []
            for k in VIQ.gene_keys:
                v = getattr(record, k)

                # gene_keys = ['rank', 'chr', 'gene', 'transcript', 'vid',
                #              'csq', 'dist', 'denovo', 'type', 'zygo', 'csn',
                #              'pldy', 'sites', 'par', 'loc', 'length', 'gqs',
                #              'gflg', 'gflpr', 'ppp', 'vpene', 'breath', 'fix',
                #              'viqscr', 'p_scor', 's_scor', 'phev_k',
                #              'vvp_svp', 'vaast', 'rprob', 'g_tag', 'p_mod',
                #              's_mod', 'g_tag_scr', 'clinvar', 'var_qual',
                #              'rid', 'loc2', 'maf', 'incndtl', 'payload']

                if k == 'pldy':
                    v = v['pldy']
                if k == 'g_tag_scr':
                    v = ','.join(v)
                if k == 'var_qual':
                    v = json.dumps(v)
                if k == 'loc2':
                    v = '{0}:{1}-{2}'.format(v['chr'], v['start'], v['end'])

                gene.append(v)
            genes.append(gene)
        
        return pd.DataFrame(data=genes, columns=VIQ.gene_keys)

#-----------------------------------------------------------------------------

def main():
    import argparse

    parser = argparse.ArgumentParser(
         description=('This python module is intended to be imported, ' +
                      'however it can be tested by passing a vIQ output ' +
                      'file, in which case it will print portions of the ' +
                      'parsed data'),
         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('file',
                        help='A vIQ output file')
    args = parser.parse_args()

    viq = VIQ(args.file)
    for record in viq.gene_records.values():
         print("\t".join([record.gene, record.viqscr, json.dumps(record.var_qual)]))

    for record in viq.mim_records.values():
         print("\t".join([record.gene, record.miqscr, record.disease]))

    for record in viq.multi_records.values():
        print("\t".join([record.phev_gene, record.mim, record.viqscr, record.miqscr]))

    ## LOH DETECTED:NO
    print('loh_detected: {}'.format(viq.meta['loh_detected']))

    ## SKEW DETECTED:NO (0)
    print('skew_detected: {}'.format(viq.meta['skew_detected']))

    ## Parental relationships MSUM:1 FSUM:1 FMR:1 DSUM:778 DSR:389
    print('parental_relationships: {}'.format(viq.meta['parental_relationships']))
    print('pr_msum: {}'.format(viq.meta['pr_msum']))
    print('pr_fsum: {}'.format(viq.meta['pr_fsum']))
    print('pr_fmr: {}'.format(viq.meta['pr_fmr']))
    print('pr_dsum: {}'.format(viq.meta['pr_dsum']))
    print('pr_dsr: {}'.format(viq.meta['pr_dsr']))

    ## SKEW DETECTED:NO (0)
    print('skew: {}'.format(viq.meta['skew']))

    ## Estimated Consanguinity RAW:5.57 % Ancestry ADJ:0 %
    print('estimated_consanguinity: {}'.format(viq.meta['estimated_consanguinity']))
    print('ec_raw: {}'.format(viq.meta['ec_raw']))
    print('ec_adj: {}'.format(viq.meta['ec_adj']))

    ## CSN PRIOR:0.5
    print('csn_prior: {}'.format(viq.meta['csn_prior']))

    ## VARIANTS_IN:14993 NUMBER OF SVs (rows):1754 PASSING_FILTERS:2307 p_obs:0.5
    print('variants_in: {}'.format(viq.meta['variants_in']))
    print('vars_snv: {}'.format(viq.meta['vars_snv']))
    print('vars_sv: {}'.format(viq.meta['vars_sv']))
    print('vars_passing: {}'.format(viq.meta['vars_passing']))
    print('vars_p_obs: {}'.format(viq.meta['vars_p_obs']))

    ## Ext. SVs in file:1657
    print('ext_svs: {}'.format(viq.meta['ext_svs']))

    ## Ext. SVs in kept:0
    print('ext_svs_kept: {}'.format(viq.meta['ext_svs_kept']))

    ## Badges trimmed:0
    print('badges_trimmed: {}'.format(viq.meta['badges_trimmed']))

    ## Badges considered:29788
    print('badges_considered: {}'.format(viq.meta['badges_considered']))

    ## Badges kept:5
    print('badges_kept: {}'.format(viq.meta['badges_kept']))

    ## Badges deleted:1
    print('badges_deleted: {}'.format(viq.meta['badges_deleted']))

    ## Proband Ancestry:Other	Relative -Log Likelihoods:	Other:10355,European(non-Finish):11019,Finish:11021,Ashkenazi:11079,Asian:11373,African:12613
    print('proband_ancestry: {}'.format(viq.meta['proband_ancestry']))
    print('pa_rlls: {}'.format(viq.meta['pa_rlls']))

    ## Proband Sex:f P(MALE):0.322425452976684
    print('proband_sex: {}'.format(viq.meta['proband_sex']))
    print('ps_prob_male: {}'.format(viq.meta['ps_prob_male']))

    ## NUM HPO TERMS:3
    print('hpo_term_count: {}'.format(viq.meta['hpo_term_count']))

    ## ADJ FOR INBREEDING:0.541174032480132
    print('adj_for_inbreeding: {}'.format(viq.meta['adj_for_inbreeding']))

    ## MODE:SINGLETON
    print('mode: {}'.format(viq.meta['mode']))

    ## Number of variants failing -e m  a:cov:120,bias:171,tot:248 x:cov:5,bias:9,tot:14
    print('num_vars_failing: {}'.format(viq.meta['num_vars_failing_em']))
    print('nvf_a: {}'.format(viq.meta['nvf_a']))
    print('nvf_x: {}'.format(viq.meta['nvf_x']))

    ## VAAST-VVP COOR:0.51
    print('vaast_vvp_coor: {}'.format(viq.meta['vaast_vvp_coor']))

    ## BLS-BND COOR:0.0049079754601227
    print('bls_bnd_coor: {}'.format(viq.meta['bls_bnd_coor']))

    ## BLS-NOA COOR:0.0049079754601227
    print('bls_nor_coor: {}'.format(viq.meta['bls_noa_coor']))

    ## PHEV-KPR COOR:0.34
    print('phev_kpr_coor: {}'.format(viq.meta['phev_kpr_coor']))

    ## CLIN-VVP-VAAST:0.725758988640042
    print('clin_vvp_vaast: {}'.format(viq.meta['clin_vvp_vaast']))

    ## COVERAGE-HOMOZYGOSITY JACC:0
    print('coverage_homozygosity_jacc: {}'.format(viq.meta['coverage_homozygosity_jacc']))

    ## COVERAGE-HETEROZYGOSITY JACC:0.04
    print('coverage_heterozygosity_jacc: {}'.format(viq.meta['coverage_heterozygosity_jacc']))

    ## PHEV-GENE-MIM COOR:0.42
    print('phev_gene_mim_coor: {}'.format(viq.meta['phev_gene_mim_coor']))

    ## K_PRIOR:0.8344
    print('k_prior: {}'.format(viq.meta['k_prior']))

    ## MIM-GEN SCORES COOR:0.08
    # print('XXX: {}'.format(viq.meta['xxx']))

    ## K_PRIOR_PROB:0.232164449818622
    # print('XXX: {}'.format(viq.meta['xxx']))

    ## U_PRIOR_PROB:0.405990586221652
    # print('XXX: {}'.format(viq.meta['xxx']))

    # SV PRIOR OBS	a	f	1	80
    print('sv_prior_obs: {}'.format(viq.meta['sv_prior_obs']))

    ## AVE DEPTH OF COVERAGE 1 MEAN:38.2462006079027 VARIANCE:266.838590703536
    print('ave_depth_of_coverage: {}'.format(viq.meta['ave_depth_of_coverage']))

    ## TYPE FREQUENCIES:	1:0.862332695984704	2:0.137667304015296	SUM:523
    print('type_frequencies: {}'.format(viq.meta['type_frequencies']))

    ## TFA (INDEL) ADJUSTMENT	STATUS:on	TFA:0.389454872134458	Pval:0.000211756350740666	exp freq:0.0878152942746563
    print('tfa_adjustment: {}'.format(viq.meta['tfa_adjustment']))
    print('tfa_status: {}'.format(viq.meta['tfa_status']))
    print('tfa_tfa: {}'.format(viq.meta['tfa_tfa']))
    print('tfa_pval: {}'.format(viq.meta['tfa_pval']))
    print('tfa_exp_freq: {}'.format(viq.meta['tfa_exp_freq']))

    ## LOW QUALITY LIST FILE DETECTED:NO. Frac Vars w/out reads:0
    print('low_quality_list_file: {}'.format(viq.meta['low_quality_list_file']))
    print('low_quality_list_file_frac: {}'.format(viq.meta['low_quality_list_file_frac']))

    ## NULL READ COUNT ADJ. With:1 Without:1
    print('null_read_count_adj_with: {}'.format(viq.meta['null_read_count_adj_with']))
    print('null_read_count_adj_without: {}'.format(viq.meta['null_read_count_adj_without']))
    # print('XXX: {}'.format(viq.meta['xxx']))

    ## CMD:/home/ubuntu/vIQ/bin/vIQ2stripe -a /home/ubuntu/fabric_viq_workflow/snakemake/rady_benchmarks/viq.config -b 0.95 -B viq_ontos/bad.txt -c  -d  -T WGS -e m -f 0.005 -g  -H  -h  -k  -l VIQ/coding_dist.CT.CS4EBE.viq_list.txt -m s -M Phevor/phv_renorm_mim.CT.CS4EBE.txt -o  -p 0.5 -Q m -q n -r n -s Phevor/phv_renorm.CT.CS4EBE.txt -v  -w  -x  -y  -z
    print('cmd: {}'.format(viq.meta['cmd']))

    ## Threshold for SV support (BFT):0.698970004336019
    print('threshold_sv_support: {}'.format(viq.meta['threshold_sv_support']))

    ## K:6.51329346539193
    print('k: {}'.format(viq.meta['k']))

    ## BEE:0
    print('bee: {}'.format(viq.meta['bee']))

    ## VERSION:7.0.02
    print('version: {}'.format(viq.meta['version']))

    ## GMT:Tue Jun 21 22:32:36 2022
    print('gmt: {}'.format(viq.meta['gmt']))

    ## EOF
    print('eof: {}'.format(viq.meta['eof']))

    # print(viq.genes_as_list())
    
if __name__ == "__main__":
    main()

