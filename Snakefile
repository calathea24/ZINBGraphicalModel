configfile: "config.yaml"

simulations=config["simulations"]
ShrinkMet=config["ShrinkMet"]
alg=config["alg"]
ps=config["ps"]
ns=range(config["ns1"], config["ns2"], config["ns3"])
degrees=config["degrees"]
graph=config["graph"]
model=config["model"]
seqDepth=config["seqDepth"]
preProcess=config["preProcess"]
pValModel=config["pValModel"]
FDRModel=config["FDRModel"]
graphSelect=config["graphSelect"]
ZI=config["ZI"]
ZImodeling=config["ZImodeling"]
iterations=range(config["iterations"])

outputs = ["results/mcc"+"_"+s+met+mod+"/"+s+"_"+prep+"_"+met+"_"+a+"_p"+str(p)+"_n"+str(n)+"_"+str(d)+"_"+g+"_"+mod+"_seqDepth"+str(dep)+"_"+pMod+"_"+FDRMod+"_"+gSe+"_ZI"+zi+"_ZImod"+zimod+"_iter"+str(i)+"_mcc.rda" for s in simulations for prep in preProcess for met in ShrinkMet for a in alg for p in ps for n in ns for d in degrees for g in graph for mod in model for dep in seqDepth for pMod in pValModel for FDRMod in FDRModel for gSe in graphSelect for zi in ZI for zimod in ZImodeling for i in iterations]

rule all:
	input: outputs

rule simulateER:
	output: temp("results/pcor/{sim}_true_{p}_{n}_{d}_{g}_{mod}_{dep}_{zi}_{i}_pcor.rda"),temp("results/simulations/{sim}_{p}_{n}_{d}_{g}_{mod}_{dep}_{zi}_{i}_data.rda")
	params:
		P = "{p}",
		N = "{n}",
		degree = "{d}",
		graph = "{g}",
		model = "{mod}",
		depth = "{dep}",
		ZI = "{zi}"
	conda: "envs/R.yaml"
	script: "snakescripts/Data_simulation.R"

rule dataPreprocessing:
	input: "results/simulations/{sim}_{p}_{n}_{d}_{g}_{mod}_{dep}_{zi}_{i}_data.rda", "results/pcor/{sim}_true_{p}_{n}_{d}_{g}_{mod}_{dep}_{zi}_{i}_pcor.rda"
	output: temp("results/preprocess/{sim}_{prep}_{p}_{n}_{d}_{g}_{mod}_{dep}_{zi}_{zimod}_{i}_data.rda")
	params:
		preProcess = "{prep}",
		model = "{mod}",
		ZI = "{zi}",
		ZImodeling = "{zimod}"
	conda: "envs/R.yaml"
	script: "snakescripts/DataPreprocessing.R"


rule pcorShrink:
	input: "results/preprocess/{sim}_{prep}_{p}_{n}_{d}_{g}_{mod}_{dep}_{zi}_{zimod}_{i}_data.rda"
	output: temp("results/pcor/{sim}_{prep}_{met}_{a}_{p}_{n}_{d}_{g}_{mod}_{dep}_{zi}_{zimod}_{i}_pcor.rda"), temp("results/pcor/{sim}_{prep}_{met}_{a}_{p}_{n}_{d}_{g}_{mod}_{dep}_{zi}_{zimod}_{i}_lambda.rda")
	params:
		ShrinkMet = "{met}",
		alg = "{a}",
		ZI = "{zi}",
		depth = "{dep}",
		ZImodeling = "{zimod}"
	conda: "envs/R.yaml"
	script: "snakescripts/Shrinkage.R"


rule pcor2edges:
	input: "results/pcor/{sim}_{prep}_{met}_{a}_{p}_{n}_{d}_{g}_{mod}_{dep}_{zi}_{zimod}_{i}_pcor.rda", "results/pcor/{sim}_{prep}_{met}_{a}_{p}_{n}_{d}_{g}_{mod}_{dep}_{zi}_{zimod}_{i}_lambda.rda"
	output: temp("results/edges/{sim}_{prep}_{met}_{a}_{p}_{n}_{d}_{g}_{mod}_{dep}_{pMod}_{FDRMod}_{gSe}_{zi}_{zimod}_{i}_edges.rda")
	params:
		ShrinkMet = "{met}",
		P = "{p}",
		N = "{n}",
		pValModel = "{pMod}",
		FDRModel = "{FDRMod}",
		graph = "{g}",
		alg = "{a}",
		model = "{mod}",
		preProcess = "{prep}",
		graphSelect = "{gSe}",
		depth = "{dep}"
	conda: "envs/R.yaml"
	script: "snakescripts/Pcor2Edges.R"

rule mcc:
	input: "results/pcor/{sim}_true_{p}_{n}_{d}_{g}_{mod}_{dep}_{zi}_{i}_pcor.rda","results/edges/{sim}_{prep}_{met}_{a}_{p}_{n}_{d}_{g}_{mod}_{dep}_{pMod}_{FDRMod}_{gSe}_{zi}_{zimod}_{i}_edges.rda"
	output: "results/mcc_{sim}{met}{mod}/{sim}_{prep}_{met}_{a}_p{p}_n{n}_{d}_{g}_{mod}_seqDepth{dep}_{pMod}_{FDRMod}_{gSe}_ZI{zi}_ZImod{zimod}_iter{i}_mcc.rda"
	params:
		P = "{p}",
		N = "{n}",
		degree = "{d}",
		graph = "{g}",
		model = "{mod}",
		ShrinkMet = "{met}",
		alg = "{a}",
		pValModel = "{pMod}",
		FDRModel = "{FDRMod}",
		i = "{i}"
	conda: "envs/R.yaml"
	script: "snakescripts/MCC.R"


#degreePlotFiles = ["results/mcc/"+s+"_"+met+"_"+str(a)+"_p"+str(p)+"_n"+str(n)+"_"+str(d)+"_"+str(g)+"_"+str(mod)+"_"+str(pMod)+"_"+str(FDRMod)+"_iter"+str(i)+"_mcc.csv" for s in simulations for met in ShrinkMet for a in alg for p in ps for n in ns for d in degrees for g in graph for mod in model for pMod in pValModel for FDRMod in FDRModel for i in iterations]

#rule makeDegreePlot:
#	input: degreePlotFiles
#	output: "results/degreePlotData.csv"
#	shell:
#		"echo \"Method,Graph,Algorithm,pVal,FDR,P,N,degree,i,MCC\" > {output};cat {input} >> {output}"

#rule degreePlot:
#	input: "results/degreePlotData.csv"
#	output: "results/plots/degreePlot.pdf"
#	script: "snakescripts/degreePlot.R"

