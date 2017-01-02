#! /usr/bin/python

#Object: Stepping stone model
#Aim: Simulating SNPs data for genom scan methods comparison
#Author: Pierre de Villemereuil
#October 2013

#************************************************************************************
#---------------------------------Header---------------------------------------------
#************************************************************************************

#Setting allele size
from simuOpt import setOptions
setOptions(alleleType='binary') #Binary for two alleles states (SNPs)

#Importing simuPOP
import simuPOP as sim
from simuPOP.utils import saveCSV
from simuPOP.sampling import drawRandomSample
from simuPOP.utils import migrSteppingStoneRates

#Importing some other python libraries
from math import *
from numpy import *
from random import *

#**********************************************************************************************
#-----------------------------Simulation parameters-------------------------------------------
#**********************************************************************************************

#Population(s) size(s)
popsize=500

#Number of generations
numgen=400
#Initial fission event
atgen=100

#Number of chromosomes
numchrom=10
#Number of loci per chrom.
numloc=500
#Create vector of loci numbers
vecloc=[numloc]*numchrom

#Migration rate (probability for one individual to disperse)
m=0.0045

#Setting up env values
vec_env=[0]*16
j=float(-1.5)
for i in range(16):
	vec_env[i]=round(j,1)
	j=j+0.2

#Vector of selected loci
locisel=[2793,1850,583,4083,3349,860,4785,706,947,939,1819,925,403,2867,2897,97,3102,2618,708,1190,2471,1533,3924,2395,2690,2926,1511,668,4826,4755,638,4148,1777,1869,2252,4326,397,3416,3171,2451,1233,2055,3013,3202,1055,3484,2984,2145,4547,4831]

#**********************************************************************************************
#-------------------------------------- Functions ---------------------------------------------
#**********************************************************************************************


#Function naming the alleles
#Depends on the number of chromosomes (chrom) and loci (loc)
def allele_naming(chrom,loc):
	res=[]
	for i in range(chrom):
		for j in range(1,loc+1):
			res.append(chr(65+i)+str(j))
	return res

#Function to attribute env values to individuals' info fields
def env_set(pop):
#Getting population size
	sim.stat(pop,popSize=True)
	subsize=pop.dvars().subPopSize
	numpop=len(subsize)
#Attribute environmental value to all individuals of the same population
	for i in range(numpop):
		pop.setIndInfo(vec_env[i],'env',subPop=i)
	return True

#Link function between env value and selection
#Modified 'logitistic' function
def fit_func(x):
# beta is the slope of the 'logistic' function
	beta=5
# coeff_s is the 'coefficient of selection'
	coeff_s=250
	res=(1-exp(-beta*x))/(1+exp(-beta*x))
	res=res/coeff_s
	return res

#Defining fitness according to selection
def fit_env(geno,env):
  #N is the number of selected loci
	N=len(geno)
  #s is the selective value (as a function of the environment)
	s=fit_func(env)
	t11=0
	t00=0
  #Calculating the number of (0,0) and (1,1) homozygotes
	for i in range(N/2):
		a1=geno[i*2]
		a2=geno[i*2+1]
		if (a1+a2==0):
			t00=t00+1
		if (a1+a2==2):
			t11=t11+1
  #w is the fitness of the individual
	w=((1+s)**t11)*((1-s)**t00)
	return w


#Function defining a constant subPop size (popsize) for any number of subpops (demographic model)
def demo(pop):
	sim.stat(pop,popSize=True)
	subsize=pop.dvars().subPopSize
#If subsize is of length 1, then it is a integer and len() does not work
	if type(subsize)==type(1):
		numpop=1
	else:
		numpop=len(subsize)
	vecsize=[popsize]*(numpop)
	return vecsize


#**********************************************************************************************
#----------------------------------------Running the simulator------------------------------
#**********************************************************************************************

num_sim=100

#Loop over the number of simulations to be done
for k in range(1,(num_sim+1)):
	print('Simulation number '+str(k))
	#Creating initial population of size popsize and "caryotype" vecloc
	pop=sim.Population(size=[popsize],loci=vecloc,infoFields=['migrate_to','fitness','env'],lociNames=allele_naming(numchrom,numloc))
	#--------------------------Main evolving process------------------------
	pop.evolve(
	#Initializing sex and genotype
		initOps=[
			sim.InitSex(),
			sim.InitGenotype(freq=[0.5,0.5]),
		],
		preOps=[
	#Fisson of population into 16 at 'atgen' generation
		sim.SplitSubPops(proportions=[0.0625]*16, at=atgen),
	#Migration using a simple stepping stone model
		sim.Migrator(rate=migrSteppingStoneRates(m,16,circular=False),begin=atgen),
	#Selection process
	#Set environmental value (env infoField) for each individual in the population
	#Takes place at each generation after the first fission
		sim.PyOperator(env_set,begin=atgen),
	#Selection occures at selected loci according to env information field
		sim.PySelector(fit_env,loci=locisel,begin=atgen),
		],
	#Mating at random (pangamy)
		matingScheme=sim.RandomMating(
	#Fixed population size (fixed at 'popsize')
		subPopSize=demo,
	#Recombination
		ops=[sim.Recombinator(rates=0.002)]
		),
		postOps=[
	#Mutation rate 10e-6
		sim.SNPMutator(u=0.000001,v=0.000001)
		],
	#Evolve for a number 'numgen' of generations
		gen = numgen
	)
	#Getting population informations (number of subpopulations, population size)
	sim.stat(pop,popSize=True)
	subsize=pop.dvars().subPopSize
	numpop=len(subsize)
	#Setting environmental value for all individuals in each subpopulation
	for i in range(numpop):
		pop.setIndInfo(vec_env[i],'env',subPop=i)
	#Sampling 20 individuals at random in each population
	sample = drawRandomSample(pop, sizes=[20]*numpop)
	#Adding population name to the field of individuals
	sample.addInfoFields('pop_name')
	vecname=[]
	for i in range(1,numpop+1):
		vecname=vecname+[i]*20
	sample.setIndInfo(vecname,'pop_name')
	#Saving the data into a .csv format
	saveCSV(sample,filename="sim"+str(k)+".csv",infoFields=['pop_name','env'],sexFormatter=None,affectionFormatter=None,header=False)

