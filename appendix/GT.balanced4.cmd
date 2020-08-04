
#NEXUS

BEGIN NETWORKS;

Network net = (((a:0.1,b:0.1):0.1,(c:0.15,d:0.15):0.05):0.05,((e:0.02,f:0.02):0.13,(g:0.04,h:0.04):0.11):0.1);


END;


BEGIN TREES;


Tree geneTree = ((((a0,((a2,a3),a1)),(((b2,b3),b1),b0)),(((c0,c1),(c2,c3)),(d0,(d1,(d2,d3))))),(((e0,(e1,(e2,e3))),(((f2,f3),f1),f0)),((((g2,g3),g1),g0),(((h2,h3),h1),h0))));


END;


BEGIN PHYLONET;

CalGTProb net (geneTree) -a <f:f2,f3,f1,f0;a:a0,a2,a3,a1;h:h2,h3,h1,h0;d:d0,d1,d2,d3;c:c0,c1,c2,c3;e:e0,e1,e2,e3;b:b2,b3,b1,b0;g:g2,g3,g1,g0> ;

END;
