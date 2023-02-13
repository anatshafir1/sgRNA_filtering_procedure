#import naiveMC as naive

class Candidate:

    #a toy example:
     #AATCCAATGTCTCTGGTTTT, 0.5, 0.19999999999999996, {'AT3G21670.1': 0.19999999999999996}, {'AT3G21670.1': [['AATCCAACGTCTCTGGTTTT', {7: ['T', 'C']}]]}

     #instead of genes_list and a targets_list: dictionaries.  in the match sites dict, use tuple insted of a set to describe the mm
    #def __init__(self, seq, fraction_of_cut = 1, cut_expectation = 1, genes_score_dict = dict(), targets_dict = dict(), lowest_cut_site_prob = 1):
    def __init__(self, seq, cut_expectation = 0, genes_score_dict = dict(), targets_dict = dict()):

        '''
        :param seq:
        #:param fraction_of_cut: the fraction of the genes in the current node, expected to be cleaved with high probability
        :param cut_expectation: the probability that ALL of the genes that expected to be cleaved by this sgRNA with high probability will be cleaved
        :param genes_score_dict: key: gene_name. value: cleaving pobability of this gene
        :param targets_dict (old name: match_sites_dict: key: gene_name. value: list of lists. each list: a match sites and a mismatches dictionary. right now it is not implemented like this at the Naive code- have to made deabuging, and to deaside maybe it is better not to save the mm locations in here.
        :param lowest_cut_site: the probability of cutting the site with the lowest cut probability by this sgRNA, among the sites of the sub tree
        width: multipication of cleaving oll of the tergets of the node in the tergets tree.
        :return:
        '''
        self.seq = seq
        #self.fraction_of_cut = fraction_of_cut
        self.cut_expectation = cut_expectation
        self.genes_score_dict = genes_score_dict
        self.targets_dict = targets_dict #change the name to missmatch site dict.
        #self.width = width
        self.score_2 = None  #the score for the objective function
        self.num_of_genes_above_thr = 0
        self.cleave_all_above_thr = 1
    '''
    def setFractionOfNodesGenes(self,genes_fract):
        self.FractionOfNodesGenes = genes_fract

    def setScore(self, score):
        self.score = score

    def set_genes_lst(self, genes_lst):
        self.genes_lst = genes_lst

    def set_targets_lst(self, targets_lst):
        self.targets_lst = targets_lst
    '''

    def fill_default_fildes(self, gene_names):
        '''
        for use when the sgRNA is constracted in the leaf of the BU tree
        :param gene_names: a list of gene names with a perfect match to the given sgRNA
        :return:
        '''
        self.genes_score_dict = dict()
        #print("gene names:" , gene_names)
        #self.fraction_of_cut = 1
        self.cut_expectation = 1
        self.lowest_cut_site_prob = 1
        #self.width = 1
        for gene_name in gene_names:
            self.genes_score_dict[gene_name] = 1
            self.targets_dict = {gene_name:[[self.seq, {}]]}


    def __str__(self):
        return self.seq  + ", " + str(self.cut_expectation) + ", " + str(self.genes_score_dict) + ", " + str(self.targets_dict)

    def __repr__(self):
        return self.__str__()
#(other_genes_score_dict, other_targets_dict[gene_name], gene_name)
    #    def add_known_site(self, other_match_sites_dict, cleaving_probability, gene_name):

    def add_known_site(self, other_targets_dict, cleaving_probability, gene_name):
        '''
        update the genes_score_dict and the match sites_dict
        :param other_targets_dict:
        :param cleaving_probability: the cleaving probability of the gene by tbe other sgRNA
        :param gene_name:
        :return:
        '''
        print("cleaving prob:", cleaving_probability)
        print("other match sites dict: ", other_targets_dict)
        if cleaving_probability < self.lowest_cut_site_prob:
            self.lowest_cut_site_prob = cleaving_probability
        #self.width = self.width * cleaving_probability
        if gene_name in self.genes_score_dict:
            ##update the cleaving probability
            self.genes_score_dict[gene_name] = 1 - (1 - self.genes_score_dict[gene_name]) * (1-cleaving_probability)
        else:
            self.genes_score_dict[gene_name] = cleaving_probability
        #now, update the match_sites_dict
        if gene_name in self.targets_dict:  #if it is in genes_score_dict, in is in here as well, but the code is flaxible like this
            self.targets_dict[gene_name] = self.targets_dict[gene_name] + other_targets_dict[gene_name]
        else:
             self.targets_dict[gene_name] = other_targets_dict[gene_name]

    def add_known_sites_of_gene(self, other_targets_dict, other_gene_cleaving_prob ,gene_name):
        '''
        :param other_targets_dict:
        :param gene_name:
        :return:
        '''
        if not gene_name in self.genes_score_dict: #still need to check the case of redundency in targets among different genes.
            self.genes_score_dict[gene_name] = other_gene_cleaving_prob
            self.targets_dict[gene_name] = other_targets_dict[gene_name]
        else:
            #print("here")
            #print(self.genes_score_dict[gene_name],"agfdaf", other_gene_cleaving_prob)
            new_cleaving_prob = 1- (1-self.genes_score_dict[gene_name])*(1-other_gene_cleaving_prob)
            self.genes_score_dict[gene_name] = new_cleaving_prob
            self.targets_dict[gene_name] = self.targets_dict[gene_name] + other_targets_dict[gene_name]


    def add_known_sites(self, other_candidate, node):
        '''
        and update the lowest_cut_site_prob
        :param other_candidate:
        :return:
        '''
        for gene_name in other_candidate.genes_score_dict.keys():
            self.add_known_sites_of_gene(other_candidate.targets_dict, other_candidate.genes_score_dict[gene_name], gene_name)



    def recalculate_cut_prob(self, Omega, len_genes_sg_dict):
        '''
        This function is sutable for the implementation with Omega. Not the current implementation
        :param Omega:
        :param len_genes_sg_dict:
        :return:
        '''
        cut_prob = 1
        num_of_cut_genes = 0
        for gene in self.genes_score_dict:
            if self.genes_score_dict[gene] > Omega:
                cut_prob = cut_prob * self.genes_score_dict[gene]
                num_of_cut_genes += 1
        self.cut_expectation = cut_prob
        self.fraction_of_cut = num_of_cut_genes / len_genes_sg_dict


    def total_num_of_mm(self):
        num_of_mm = 0
        for val in self.targets_dict.values(): #here the value is a list of targets. Each represented by list. Each list: a match sites and a mismatches dictionary.
            for target in val:
                #for mm in target[1].values():
                num_of_mm += len(target[1])
        return  num_of_mm
