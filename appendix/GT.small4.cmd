
#NEXUS

BEGIN NETWORKS;

Network net = (((a:1.0,b:1.0):0.5,c:1.5):0.3,(d:0.9,e:0.9):0.9);


END;


BEGIN TREES;


Tree geneTree = ((((((a2,a3),a1),a0),(b0,((b2,b3),b1))),((c2,c3),(c0,c1))),(((d2,d3),(d0,d1)),(e0,(e1,(e2,e3)))));


END;


BEGIN PHYLONET;

CalGTProb net (geneTree) -a <c:c2,c3,c0,c1;e:e0,e1,e2,e3;a:a2,a3,a1,a0;b:b0,b2,b3,b1;d:d2,d3,d0,d1> ;

END;
