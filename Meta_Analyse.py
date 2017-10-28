#!/usr/bin/env python
#-*- coding: utf-8 -*-

import pandas
import numpy as np
from math import exp, log, sqrt
from scipy.stats import chi2

""" Jeux de données pour résultats d'études disponibles en nombre
d'événements par essai  """
data_ev = pandas.read_table("Data_Ev.txt",sep = '\t',header = 0)


""" Jeux de données pour résultats disponibles directement en risque relatif
effect (Risque relatif) """
data_TE = pandas.read_table("Data_TE.txt",sep = '\t',header = 0)


class Meta:
	"""
	classe meta-analyse regroupe toutes les fonctions:
	- d'entrée des données
	- de calcul
	- de représentation graphique
	Les calculs de la méta analyse se bases sur la méthode du logarithme du risque relatif.
	"""
	
	def __init__(self, studLab=None, x1=None, n1=None, x0=None, n0=None, RR=None, bi=None, bs=None, 
	data_type=None):
			self.studLab = studLab
			self.x1 = x1
			self.n1 = n1
			self.x0 = x0
			self.n0 = n0
			self.RR = RR
			self.bi = bi
			self.bs = bs
			self.data_type = data_type
		
	@classmethod	
	def events(cls, studLab, x1, n1, x0, n0):
		""" Entrée des données pour les résultats d'études disponibles en nombre
		d'événements par essai.
		Les paramètres d'entrée sont :
			- studLab : Nom de l'essai
			- x1 : Nombre d'événements  du groupe traité 
			- n1 : Effectif du groupe traité
			- x0 : Nombre d'événements du groupe contrôle
			- n0 : Effectif du groupe contrôle
		"""
		myMeta =  Meta(studLab=studLab, x1=x1, n1=n1, x0=x0, n0=n0, data_type = 1)
		return (myMeta)
		
			
	@classmethod	
	def TE(cls, studLab, RR, bi, bs):
		""" Entrée des données pour les résultats disponibles directement en risque relatif 
		Les paramètres d'entrée sont :
			- studLab : Nom de l'essai
			- RR : Risque relatif
			- bi : Borne inférieur de l'interval de confiance à 95 % du risque relatif 
			- bs : Borne supérieur de l'interval de confiance à 95 % du risque relatif 
		"""
		myMeta =  Meta(studLab = studLab , RR = RR, bi = bi, bs = bs, data_type = 2 )
		return (myMeta)
	

	def Heterogeneity_Test(self):
		""" Cette Méthode effectue le test d'hétérogénéité de la méta-analyse.
		Elle retourne dans le dictionnaire "self.heterogeneity_test" :
			- Q : Statistique du test
			- ddl : Le nombre de degré de liberté du test
			- pvalue : La pvalue du test
		"""
		
		# Statistique du Test
		Q = sum((self.logRR)**2 *self.w) - (sum(self.logRR*self.w)**2/sum(self.w))
		# Nombre de ddl du test
		df = len(self.studLab) - 1
		# P value du test
		pvalue = 1 - chi2.cdf(Q, df)
		self.heterogeneity_test = {'Q' : round(Q,3), 'df' : df,'pvalue': round(pvalue,6)}


	def Assocation_Test(self):
		""" Cette Méthode effectue le test d'association de la méta-analyse.
		Elle retourne dans le dictionnaire "self.association_test" :
			- Z : Statistique du test
			- pvalue : La pvalue du test
		"""
		Z =  sum(self.logRR*self.w)**2/sum(self.w)
		pvalue = 1 - chi2.cdf(Z, 1)	
		self.association_test = {'Z' : round(Z,3), 'pvalue': round(pvalue,6)}
		
		
		
	def Commun_TE(self):
		""" Cette Méthode calcule l'effet du traitement commun.
		Elle retourne dans le dictionnaire "self.commun_TE" :
			- 'Common Treatement Effect ' : L'estimation ponctuelle de  l'effet du traitement commun
			- Borne Inf' : La borne inférieur de l'interval de confiance
			- Borne Sup' : La borne supérieur de l'interval de confiance
		"""
		
		TE_Com = exp(sum(self.logRR*self.w)/sum(self.w))
		var_log_theta_com = 1/sum(self.w)
		bi =  exp(log(TE_Com) - 1.96*sqrt(1/sum(self.w)))
		bs =  exp(log(TE_Com) +1.96*sqrt(1/sum(self.w)))
		self.commun_TE = {'Common Treatement Effect ' : round(TE_Com,3), 'Borne Inf':
		round(bi,3), 'Borne Sup': round(bs,3)}
		




	def MetaAnalyse(self):

		""" Méthode permettant d'effectuer l'ensemble des calculs de la méta-analyse:
			- Le test d'hétérogénéité
			- L'estimation de l'effet traitement commun
			- Le test d'association
			
		
		Etape 1: Estimation des variables intermédiaires permettant les calculs:
			- logRR : Logarithme du risque relatif de chaque essai (Effet du Traitement) .
			- varLogRR : Variance de logRR .
			- w: Inverse de la variance de logRR (Estimateur du poids de l'essai
			dans le calcul de l'effet traitement commun . """
	

		""" Calcul de logRR et varLogRR  à partir des  données directement en risque relatif  """
		if  self.data_type == 2:
			self.logRR = np.log(self.RR)
			self.varLogRR = ((np.log(self.bs) - np.log(self.bi))/2/1.96)**2
			
			
		#Calcul de logRR et varLogRR  à partir des données disponibles en nombre
		#d'événements par essai 
		else:
			self.logRR = np.log((self.x1/self.n1)/(self.x0/self.n0))
			self.varLogRR = 1/self.x1 -  1/self.n1 +  1/self.x0 -  1/self.n0
			
		self.w = 1/self.varLogRR

		"""Etape 2 : Calculs """
		self.Heterogeneity_Test()
		self.Commun_TE()
		self.Assocation_Test()



""" Fonction Principale  """
if __name__ == '__main__':
	
	""" I) Entrée des données pour les 2 versions """
	myMeta = Meta.events(data_ev.Essai, data_ev.x1, data_ev.n1, data_ev.x0, data_ev.n0)
	myMeta = Meta.TE(data_TE.Essai, data_TE.RR, data_TE.bi, data_TE.bs)
	
	""" Calculs de la méta-analyse"""
	myMeta.MetaAnalyse()

	print(myMeta.commun_TE)
	print(myMeta.heterogeneity_test)
	print(myMeta.association_test)
	







