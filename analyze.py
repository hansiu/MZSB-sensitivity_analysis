import numpy as np
import matplotlib.pylab as pl
import pickle
import sys


def Morris():
    N = 1000
    an_m = 'morris'
    smp_m = 'morris'

    sa_results = pickle.load(open("SA_results_" + str(N) + "_" + an_m + "_" + smp_m, "rb"))

    data_mu = []
    data_mu_star = []
    data_sigma = []

    banned = ['default_size', 'c1_size', 'cell_size', '_c1_iron_in_Plasma_0']

    for out in sa_results.keys():
        mu_row = []
        mu_star_row = []
        sigma_row = []
        for i, n in enumerate(sa_results[out]['names']):
            if n not in banned:
                mu_row.append(sa_results[out]['mu'][i])
                mu_star_row.append(sa_results[out]['mu_star'][i])
                sigma_row.append(sa_results[out]['sigma'][i])
        data_mu.append(mu_row)
        data_mu_star.append(mu_star_row)
        data_sigma.append(sigma_row)

    columns = [x.split('_k1')[0] for x in sa_results['iron_in_Plasma_c1']['names'] if x not in banned]
    rows = [k.split('_')[2] for k in sa_results.keys()]

    n_rows = len(data_mu)
    index = np.arange(len(columns))
    bar_width = 0.8
    for data in [data_mu_star, data_sigma]:
        colors = pl.cm.tab20(np.linspace(0, 1, n_rows))
        y_offset = np.zeros(len(columns))
        pl.figure(figsize=(30, 8))
        cell_text = []
        for row in range(n_rows):
            pl.bar(index, data[row], bar_width, bottom=y_offset, color=colors[row])
            y_offset = y_offset + data[row]
            cell_text.append(['%.3f' % x for x in data[row]])
        # Reverse colors and text labels to display the last value at the top.
        colors = colors[::-1]
        cell_text.reverse()
        the_table = pl.table(cellText=cell_text,
                             rowLabels=rows,
                             rowColours=colors,
                             colLabels=columns,
                             loc='bottom')
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(9)
        pl.subplots_adjust(bottom=0.2)
        pl.xticks([])
        pl.margins(x=0)
        if data == data_mu_star:
            pl.savefig('Morris_mu_star_results.png', bbox_inches="tight", dpi=300)
        elif data == data_sigma:
            pl.savefig('Morris_sigma_results.png', bbox_inches="tight", dpi=300)
        pl.clf()


def Sobol():
    N = 1000
    an_m = 'sobol'
    smp_m = 'saltelli'

    sa_results = pickle.load(open("SA_results_" + str(N) + "_" + an_m + "_" + smp_m, "rb"))
    names = ['r1_k1', 'r19_k1', 'r11_k1', 'r6_k1', 'r23_k1', 'r12_k1', 'r24_k1', 'r8_k1', 'r17_k1', 'r10_k1', 'r26_k1',
             'r15_k1', 'r21_k1', 'r20_k1', 'r7_k1', 'r29_k1', 'r25_k1', 'r9_k1', 'r18_k1', 'r27_k1', 'r16_k1', 'r22_k1',
             'r4_k1', 'r14_k1', 'r13_k1', 'r28_k1', 'r2_k1', 'r5_k1', 'r3_k1', 'default_size', 'c1_size', 'cell_size',
             '_c1_iron_in_Plasma_0']

    data_s1 = []
    data_st = []
    data_s2 = []

    banned = ['default_size', 'c1_size', 'cell_size', '_c1_iron_in_Plasma_0']

    for out in sa_results.keys():
        s1_row = []
        st_row = []
        s2_row = []
        for i, n in enumerate(names):
            if n not in banned:
                s1_row.append(sa_results[out]['S1'][i])
                st_row.append(sa_results[out]['ST'][i])
                s2_row.append(sa_results[out]['S2'][i, :29])

        data_s1.append(s1_row)
        data_st.append(st_row)
        data_s2.append(s2_row)

    columns = [x.split('_k1')[0] for x in names if x not in banned]
    rows = [k.split('_')[2] for k in sa_results.keys()]

    n_rows = len(data_s1)
    index = np.arange(len(columns))
    bar_width = 0.8
    for data in [data_s1, data_st]:
        colors = pl.cm.tab20(np.linspace(0, 1, n_rows))
        y_offset = np.zeros(len(columns))
        pl.figure(figsize=(30, 8))
        cell_text = []
        for row in range(n_rows):
            pl.bar(index, data[row], bar_width, bottom=y_offset, color=colors[row])
            y_offset = y_offset + data[row]
            cell_text.append(['%.3f' % x for x in data[row]])
        # Reverse colors and text labels to display the last value at the top.
        colors = colors[::-1]
        cell_text.reverse()
        the_table = pl.table(cellText=cell_text,
                             rowLabels=rows,
                             rowColours=colors,
                             colLabels=columns,
                             loc='bottom')
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(9)
        pl.subplots_adjust(bottom=0.2)
        pl.xticks([])
        pl.margins(x=0)
        if data == data_s1:
            pl.savefig('Sobol_S1_results.png', bbox_inches="tight", dpi=300)
        elif data == data_st:
            pl.savefig('Sobol_ST_results.png', bbox_inches="tight", dpi=300)
        pl.clf()

    for i in range(n_rows):
        row_1 = data_s1[i]
        row_T = data_st[i]
        for n in range(len(row_1)):
            if row_T[n] > row_1[n]:  # possible higher order interactions
                interactions = list(np.nonzero(np.nan_to_num(data_s2[i][n]) > 0.05)[0])
                if interactions:
                    print(rows[i], columns[n], row_T[n], row_1[n], list(np.array(columns)[interactions]),
                          list(np.array(data_s2[i][n])[interactions]))


Morris()
Sobol()
