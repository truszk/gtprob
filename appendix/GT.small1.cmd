
#NEXUS

BEGIN NETWORKS;

Network net = (((a:1.0,b:1.0):0.5,c:1.5):0.3,(d:0.9,e:0.9):0.9);


END;


BEGIN TREES;


Tree geneTree = (((a0,b0),c0),(d0,e0));


END;


BEGIN PHYLONET;

CalGTProb net (geneTree) -a <c:c0;e:e0;a:a0;b:b0;d:d0> ;

END;
