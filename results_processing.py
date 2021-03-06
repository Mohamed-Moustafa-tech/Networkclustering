#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import mygene
import networkx as nx
import numpy as np
import seaborn as sns;
from numpy.core._multiarray_umath import ndarray

sns.set(color_codes=True)
import gseapy


class results_analysis():
    '''
        Performs analysis over the output of BiGAnts algorithm

        Attributes:
        -----------
        solution - the output file of BiGAnts.run_search() function
        labels - data preprocessing labels from data_preprocessing() function
        convert - indicates if gene IDs should be converted to gene names
        for the further results analysis (default - False)
        origID - indicates the original gene ids used. This field is mandatory for the enrichment analysis.
        Possible values:
            'entrezgene', 'ensembl.gene', 'symbol', 'refseq', 'unigene', etc
            for all possibe option please check  the reference for MyGene.info gene query web service
            http://docs.mygene.info/en/latest/doc/query_service.html#available_fields
    '''

    def __init__(self, nodes, labels_ids, labels, n, convert=False, origID=None):
        self.labels_ids = labels_ids
        self.nodes = nodes
        self.genes = [str(self.labels_ids[x]) for x in nodes]
        patients = np.arange(len(labels))
        patients = patients + n
        self.patients1 = [labels_ids[patients[i]] for i in range(len(labels)) if labels[i] == 1]
        self.patients2 = [labels_ids[patients[i]] for i in range(len(labels)) if labels[i] == 0]
        self.convert = convert
        self.origID = origID
        self.pts1 = [patients[i] for i in range(len(labels)) if labels[i] == 1]
        self.pts2 = [patients[i] for i in range(len(labels)) if labels[i] == 0]

        if convert == True:
            assert origID != None, "Please specify the original gene ID or set 'convert' to False"
            all_genes = self.genes
            mg = mygene.MyGeneInfo()
            out = mg.querymany(all_genes, scopes=self.origID, fields='symbol', species='human', verbose=False)
            mapping = dict()
            rev_mapping = dict()
            for line in out:
                try:
                    rev_mapping[line["symbol"]] = line["query"]
                    mapping[line["query"]] = line["symbol"]
                except KeyError:
                    print("{0} was not mapped to any gene name".format(line["query"]))
                    mapping[line["query"]] = line["query"]
                    rev_mapping[line["query"]] = line["query"]

            self.mapping = mapping

    def show_networks(self, GE, G, output=None):
        '''
        Shows the resulting subnetworks coloured wrt to their difference in expression patterns in patients subgroups

        Attributes:
        -----------
        GE - processed gene expression data from data_preprocessing() function
        G - processed PPI network from data_preprocessing() function
        output - str or PathLike or file-like object (png, eps, pdf, etc)
        '''
        # relabel solution IDs to the actual IDs

        all_genes_entr = self.genes
        all_genes = self.nodes
        if self.convert:
            all_genes_names = [self.mapping[x] for x in self.genes]

        # relabel expression matrix and the graph to the actual patients ids and gene names
        G_small = nx.subgraph(G, all_genes)
        G_small = nx.relabel_nodes(G_small, self.labels_ids)
        GE_small = GE.loc[all_genes]
        if self.convert:
            GE_small.index = all_genes_names
            G_small = nx.relabel_nodes(G_small, self.mapping)
        else:
            GE_small.index = all_genes_entr

        # compute difference in expression for each gene in different patients groups
        if self.convert:

            means1: ndarray = np.mean(GE_small[self.pts1].loc[all_genes_names], axis = 1)
            means2: ndarray = np.mean(GE_small[self.pts2].loc[all_genes_names], axis = 1)

        else:
            means1 = np.mean(GE_small[self.pts1].loc[all_genes_entr], axis = 1)
            means2 = np.mean(GE_small[self.pts2].loc[all_genes_entr], axis = 1)
        if np.mean(means1) > np.mean(means2):
            means = means1 - means2
        else:
            means = means2 - means1
        # set plotting settings
        plt.rc('font', size=20)  # controls default text sizes
        plt.rc('axes', titlesize=20)  # fontsize of the axes title
        plt.rc('axes', labelsize=20)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize=15)  # fontsize of the tick labels
        plt.rc('ytick', labelsize=25)  # fontsize of the tick labels
        plt.rc('legend', fontsize=30)
        fig = plt.figure(figsize=(15, 15))
        #        vmin = -2
        #        vmax = 2
        #        cmap = plt.cm.coolwarm(np.linspace(-0.45,0.8,20))
        #        cmap = mpl.colors.ListedColormap(cmap[10:,:-1])
        pos = nx.spring_layout(G_small)
        nx.draw_networkx_edges(G_small, pos)
        nc1 = nx.draw_networkx_nodes(G_small, pos=pos, node_color=means, node_size=1700, alpha=.7)
        nx.draw_networkx_labels(G_small, pos, font_size=22, font_weight="heavy")
        plt.colorbar(nc1)
        plt.axis('off')
        fig.tight_layout()
        # save if required
        if output != None:
            plt.savefig(output, dpi=300)
        plt.show()
    def cor_map(self, GE, output = None):
        all_genes = self.nodes
        if self.convert:
            all_genes_names = [self.mapping[x] for x in self.genes]
        else:
            all_genes_names = self.genes

        # relabel expression matrix and the graph to the actual patients ids and gene names
        GE_small = GE[self.pts1 + self.pts2].loc[all_genes]
        GE_small.index = all_genes_names
        part1 = np.corrcoef(GE_small[self.pts1])
        m1 = np.mean(np.abs(part1))
        part2 = np.corrcoef(GE_small[self.pts2])
        m2 = np.mean(np.abs(part2))
        plt.figure(figsize=(10, 5))
        plt.subplot(1, 2, 1)
        plt.title("Partition 1 mean absolute correlation {0}".format(np.round(m1,2)))
        sns.heatmap(pd.DataFrame(part1, index = all_genes_names, columns = all_genes_names), vmin = -1, vmax = 1)
        plt.subplot(1, 2, 2)
        plt.title("Partition 2 mean absolute correlation {0}".format(np.round(m2,2)))
        sns.heatmap(pd.DataFrame(part2, index = all_genes_names, columns = all_genes_names), vmin = -1, vmax = 1)
        if output != None:
            plt.savefig(output, dpi=300)
        plt.show()



    def show_clustermap(self, GE, G, true_labels=None, output=None, class_names = []):
        '''
        Shows a clustermap of the achieved solution alone or also along with the known patients groups

        Attributes:
        -----------
        GE - processed gene expression data from data_preprocessing() function
        G - processed PPI network from data_preprocessing() function
        output - str or PathLike or file-like object (png, eps, pdf, etc)
        true_labels
        '''

        #        if true_labels !=None:
        #            patients = self.patients1 + self.patients2
        #            true_patients = true_labels[0] + true_labels[1]
        #            if len(set(patients).difference(set(true_patients))) != 0:
        #                print("WARNING: Patients ids in true_labels do not match, comparisson wil not be performed")
        #                true_labels = None

        all_genes = self.nodes
        if self.convert:
            all_genes_names = [self.mapping[x] for x in self.genes]
        else:
            all_genes_names = self.genes

        # relabel expression matrix and the graph to the actual patients ids and gene names
        GE_small = GE[self.pts1 + self.pts2].loc[all_genes]
        GE_small.index = all_genes_names
        # prepare the clustermap
        p_num = GE_small.columns
        if true_labels != None:

            grouping_p_true = []
            patients1_true = true_labels[0]
            patients2_true = true_labels[1]

            for p in p_num:
                if self.labels_ids[p] in patients1_true:
                    grouping_p_true.append(5)
                elif self.labels_ids[p] in patients2_true:
                    grouping_p_true.append(6)
            grouping_p_true = pd.DataFrame(grouping_p_true, index=p_num)
            grouping_p_true.columns = ["true"]
            species = grouping_p_true["true"]
            lut = {5: '#F3FF33', 6: 'm'}
            row_colors2 = species.map(lut)

        grouping_p_clust = []
        for p in p_num:
            if self.labels_ids[p] in self.patients1:
                grouping_p_clust.append(1)
            else:
                grouping_p_clust.append(2)
        grouping_p_clust = pd.DataFrame(grouping_p_clust, index=p_num)
        grouping_p_clust.columns = ["clusters"]
        species = grouping_p_clust["clusters"]
        lut = {1: '#4FB6D3', 2: '#22863E'}
        row_colors1 = species.map(lut)

        plt.rc('font', size=5)  # controls default text sizes
        plt.rc('axes', titlesize=20)  # fontsize of the axes title
        plt.rc('axes', labelsize=20)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize=20)  # fontsize of the tick labels
        plt.rc('ytick', labelsize=20)  # fontsize of the tick labels
        plt.rc('legend', fontsize=20)

        if true_labels != None:
            g = sns.clustermap(GE_small.T, row_colors=[row_colors1, row_colors2], row_cluster=False, figsize=(15, 10),
                               cbar_kws=dict(ticks=[5, 0, -5]),
                               cmap="Spectral", yticklabels=False)

            if len(class_names) != 2:
                values = ["true class1", "true class2", "cluster1", "cluster2"]
            else:
                values = [class_names[0], class_names[1], "cluster1", "cluster2"]

            colors = ['#F3FF33', 'm', '#4FB6D3', '#22863E']
            for i in range(len(values)):
                l = values[i]
                c = colors[i]
                g.ax_col_dendrogram.bar(0, 0, color=c,
                                        label=l, linewidth=0)
            g.ax_col_dendrogram.legend(loc="upper center", ncol=2, bbox_to_anchor=(0.72, 0.87),
                                       borderaxespad=0.)
        else:
            g = sns.clustermap(GE_small.T, row_colors=[row_colors1, row_colors2], row_cluster=False, figsize=(15, 10),
                               cbar_kws=dict(ticks=[5, 0, -5]),
                               cmap="Spectral", yticklabels=False)
            g.ax_row_dendrogram.set_visible(False)
            g.cax.set_visible(False)

            values = ["cluster1", "cluster2"]
            colors = ['#4FB6D3', '#22863E']
            for i in range(len(values)):
                l = values[i]
                c = colors[i]
                g.ax_col_dendrogram.bar(0, 0, color=c,
                                        label=l, linewidth=0)
            g.ax_col_dendrogram.legend(loc="upper center", ncol=2, bbox_to_anchor=(0.72, 0.87),
                                       borderaxespad=0.)

        ax = g.ax_heatmap
        ax.set_xlabel("Genes")
        ax.set_ylabel("Patients")
        if output != None:
            g.savefig(output, dpi=300)
        plt.show()

    def jaccard_index(self, true_labels) -> object:
        """

        :rtype: object
        """

        def jac(x, y):
            if len(x) > 0 and len(y) > 0:
                return len(set(x).intersection(set(y))) / len((set(x).union(set(y))))
            else:
                return (0)

        def jac_matrix(true, pred):
            res = np.zeros((len(true), len(true)))
            for i in range(len(true)):
                for j in range(len(true)):
                    res[i, j] = jac(true[i], pred[j])
            cand1 = (res[0][0], res[1][1])
            cand2 = (res[0][1], res[1][0])
            if sum(cand1) > sum(cand2):
                return (cand1)
            else:
                return (cand2)

        ids = jac_matrix([self.patients1, self.patients2], true_labels)

        return(ids)

    def enrichment_analysis(self, library, output):
        '''
        Saves the results of enrichment analysis

        Attributes:
        -----------
        library - Enrichr library to be used. Recommendations:
            - 'GO_Molecular_Function_2018'
            - 'GO_Biological_Process_2018'
            - 'GO_Cellular_Component_2018'
            for more options check available libraries by typing gseapy.get_library_name()

        output - directory name where results should be saved
        '''
        libs = gseapy.get_library_name()
#        assert library in libs, "the library is not available, check gseapy.get_library_name() for available options"
        assert (self.convert == True) or (
                self.origID == "symbol"), "EnrichR accepts only gene names as an input, thus please set 'convert' to True and indicate the original gene ID"

        all_genes_names = [self.mapping[x] for x in self.genes]
        print(all_genes_names)
        res = gseapy.enrichr(gene_list=all_genes_names, description='pathway', gene_sets=library, cutoff=0.05, outdir=output)
        return(res.results)

    def convergence_plot(self, scores, output=None):
        '''
        Shows the convergence plot

        Attributes:
        -----------
        scores - the output of run_search() function

        output - directory name where results should be saved
        '''

        count_big, scores, avs = scores
        plt.figure(figsize=(10, 6))

        sns.set(style="whitegrid")
        plt.rc('font', size=13)  # controls default text sizes
        plt.rc('axes', titlesize=13)  # fontsize of the axes title
        plt.rc('xtick', labelsize=13)  # fontsize of the tick labels
        plt.rc('ytick', labelsize=13)  # fontsize of the tick labels
        plt.rc('legend', fontsize=13)
        zippedList = list(zip(scores, avs))
        wg = pd.DataFrame(zippedList, columns=["best score", "average score"])
        ax = sns.lineplot(data=wg, palette="tab10", linewidth=2.5)
        ax.set(xlabel="Iterations")
        ax.set(ylabel="Score")

        if output != None:
            plt.savefig(output)
        plt.show()

