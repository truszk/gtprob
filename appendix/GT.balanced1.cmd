
#NEXUS

BEGIN NETWORKS;

Network net = (((a:0.1,b:0.1):0.1,(c:0.15,d:0.15):0.05):0.05,((e:0.02,f:0.02):0.13,(g:0.04,h:0.04):0.11):0.1);


END;


BEGIN TREES;


Tree geneTree = (((a0,b0),(c0,d0)),((e0,f0),(g0,h0)));


END;


BEGIN PHYLONET;

CalGTProb net (geneTree) -a <f:f0;a:a0;h:h0;d:d0;c:c0;e:e0;b:b0;g:g0> ;

END;
