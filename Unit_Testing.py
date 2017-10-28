#tests unitaires

from Meta_Analyse import *

#test unitaire pour Méta-analyses pour critères binaires

#cas simple
myMeta = Meta.events(data_ev.Essai, data_ev.x1, data_ev.n1, data_ev.x0, data_ev.n0)
myMeta = Meta.TE(data_TE.Essai, data_TE.RR, data_TE.bi, data_TE.bs)
myMeta.MetaAnalyse()

assert myMeta.commun_TE == {'Borne Inf': 0.651, 'Borne Sup': 1.016, 'Common Treatement Effect ': 0.813}
assert myMeta.heterogeneity_test == {'Q': 0.726, 'df': 5, 'pvalue': 0.981506}
assert myMeta.association_test == {'Z': 3.311, 'pvalue': 0.06881}
