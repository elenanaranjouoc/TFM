import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, r

r('library(dplyr)')
r('library(Biobase)')
r('library(genefilter)')
r('library(hgu133plus2.db)')
r('library(multtest)')
r('library(Rgraphviz)')
r('library(GOstats)')
r('library(Category)')
r('library(stringr)')  
r('library(GSEABase)')  
r('library(GSA)')  


def grupo_genes(expresionSet):
    """
    Using the rpy2 library, we access the R terminal to create sets of genes and use the fisher test to obtain which are related to the chosen phenotypic variable.

    Parameters:
    expresionSet: R object of class expressionSet 
    """
     
    pandas2ri.activate()

    ro.globalenv["expresionSet"] = expresionSet

    r("""
    conj_genes <- function(Set, colu) {
        y = pData(Set)[,colu]
        y <-factor(y)
        if (length(unique(y)) == 2) {
            tt = rowttests(Set,y)
        }
        else {tt = rowFtests(Set,y)}
        p0 = tt$p.value
        p1 = mt.rawp2adjp(p0, "BH")
        orden.original = order(p1$index)
        p.BH = p1$adjp[orden.original, 2]
        significativos = which(p.BH < 0.05)

        G1.entreizd = unlist(mget(featureNames(Set), hgu133plus2ENTREZID))
        if (anyDuplicated(G1.entreizd) != 0) {
            set.iqr = apply(exprs(Set), 1, IQR)
            uniqGenes = findLargest(featureNames(Set), set.iqr, "hgu133plus2")
            Set1 = Set[uniqGenes, ]
            G2.entrezid = unlist(mget(featureNames(Set), hgu133plus2ENTREZID))
            G1.entrezid = G2.entrezid
        }

        seleccionados = unlist(mget(featureNames(Set[significativos, ]), hgu133plus2ENTREZID))


        # TEST FISHER
        params = new("GOHyperGParams", geneIds = seleccionados, universeGeneIds= G1.entreizd, annotation = "hgu133plus2", ontology = "BP", pvalueCutoff = 0.001, conditional = FALSE, testDirection = "over")
        overRepresented = hyperGTest(params)

        # Resultados en html
        old_wd <- getwd()
        setwd(paste0(getwd(),"/results/grupo_genes"))

        htmlReport(overRepresented, file = paste(str_remove_all(colu," "), "_overRepresented.html", sep=""))
        
        #Grafo

        png(file=paste(str_remove_all(colu," "), "_overRepresented.png", sep=""))
        plot(goDag(overRepresented))
        dev.off()
        setwd(old_wd)

    }
    conj_genes(expresionSet, "ER")
    conj_genes(expresionSet, "HER2")
    conj_genes(expresionSet, "node status")
    conj_genes(expresionSet, "tumor type")
    conj_genes(expresionSet, "B-R grade")
    """)


def gsa(expresionSet):
    """
    Using the rpy2 library, we access the R terminal to create sets of genes and use the fisher test to obtain which are related to the chosen phenotypic variable via GSA.

    Parameters:
    expresionSet: R object of class expressionSet 
    """
     
    pandas2ri.activate()

    ro.globalenv["expresionSet"] = expresionSet

    r("""
    obtener_grupos_tamaño <- function(Set, num_genes) {
        gsc = GeneSetCollection(Set, setType=GOCollection())
        names(gsc) = unlist(lapply(gsc, setName))

        old_wd <- getwd()
        setwd(paste0(getwd(),"/results/gsa"))
        save(gsc, file = "gsc.rda")
        setwd(old_wd)

        # Grupos con más de num_genes genes
        grupos_mas = gsc[which(sapply(geneIds(gsc),length) > num_genes)]
        return (grupos_mas)
    }
    gsa_f <- function(Set, gsc, colu) {
        num_node = factor(pData(Set)[, colu])
        num_node <-as.numeric(num_node)

        if (length(unique(num_node)) == 2) {
            resp = "Two class unpaired"
        }
        else {resp = "Multiclass"}

        gsa.n = GSA(exprs(Set), num_node, genenames=featureNames(Set), genesets = geneIds(gsc), resp.type=resp, nperms=1000)
        old_wd <- getwd()
        setwd(paste0(getwd(),"/results/gsa"))

        png(file=paste(colu, "_gsa.png", sep=""))
        GSA.plot(gsa.n)
        dev.off()
        

        #Asociación negativa
        as_neg = which(gsa.n$pvalues.lo < 0.05)
        neg=names(gsc[as_neg])
        write.csv(neg, file = paste("neg_",colu,".csv", sep=""))
        

        #Asociación positivaS
        as_pos = which(gsa.n$pvalues.hi < 0.05)
        pos=names(gsc[as_pos])
        write.csv(pos, file = paste("pos_",colu,".csv", sep=""))

        setwd(old_wd)

    }
    grupos <- obtener_grupos_tamaño(expresionSet, 50)
    gsa_f(expresionSet, grupos, "ER")
    gsa_f(expresionSet, grupos, "HER2")
    gsa_f(expresionSet, grupos, "node status")
    gsa_f(expresionSet, grupos, "B-R grade")
    gsa_f(expresionSet, grupos, "tumor type")
    """)
