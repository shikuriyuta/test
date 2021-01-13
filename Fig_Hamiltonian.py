import daft

#para
H = 8
V = 4.2
w1 = 0.05
w2 = 0.5
h1 = w1 + (H - w1) / 6 * 0.5
h2 = w1 + (H - w1) / 6 * 1.5
h3 = w1 + (H - w1) / 6 * 2.5
h4 = w1 + (H - w1) / 6 * 3.5
h5 = w1 + (H - w1) / 6 * 4.5
h6 = w1 + (H - w1) / 6 * 5.5
v0 = 4
v1 = 3
v2 = 2
v3 = 1
r = 0.4
r_ = 0.7

# Object
pgm = daft.PGM(shape=[H, V], node_unit = 1)

# Plate 
pgm.add_plate(daft.Plate([0, v1-0.5*r_, H-0.2, v2-v3+r_], label_offset =[5,5], label=r"$cycle$", fontsize=15))
pgm.add_plate(daft.Plate([0, v3-0.5*r_, H-0.2, v2-v3+r_], label_offset =[5,5], label=r"$score$", fontsize=15))

# Nodes
pgm.add_node(daft.Node("r12", r"$r_{1,2}$", h1, v0))
pgm.add_node(daft.Node("r13", r"$r_{1,3}$", h2, v1))
pgm.add_node(daft.Node("r14", r"$r_{1,4}$", h3, (v0+v1)/2+0.15))
pgm.add_node(daft.Node("r23", r"$r_{2,3}$", h4, v1))
pgm.add_node(daft.Node("r24", r"$r_{2,4}$", h5, v0))
pgm.add_node(daft.Node("r34", r"$r_{3,4}$", h6, (v0+v1)/2+0.15))
pgm.add_node(daft.Node("p11", r"$p_{1,1}$", w2+h1, v2))
pgm.add_node(daft.Node("p21", r"$p_{2,1}$", w2+h1, v3))
pgm.add_node(daft.Node("q11", r"$q_{1,1}$", w2+h2, (v2+v3)/2))
pgm.add_node(daft.Node("p12", r"$p_{1,2}$", w2+h3, (v2+v3)/2))
pgm.add_node(daft.Node("p13", r"$p_{1,3}$", w2+h4, (v2+v3)/2))
pgm.add_node(daft.Node("q13", r"$q_{1,3}$", w2+h5, (v2+v3)/2))

# edge
pgm.add_edge("p11", "r12", directed=False)
pgm.add_edge("p11", "r13", directed=False)
pgm.add_edge("p21", "r13", directed=False)
pgm.add_edge("q11", "r14", directed=False)
pgm.add_edge("p12", "r23", directed=False)
pgm.add_edge("p13", "r13", directed=False)
pgm.add_edge("q13", "r34", directed=False)
pgm.add_edge("p11", "p21", directed=False)
pgm.add_edge("p11", "q11", directed=False)
pgm.add_edge("p21", "q11", directed=False)
pgm.add_edge("p13", "q13", directed=False)
pgm.add_edge("r12", "r13", directed=False)
pgm.add_edge("r12", "r14", directed=False)
pgm.add_edge("r12", "r23", directed=False)
pgm.add_edge("r12", "r24", directed=False)
pgm.add_edge("r13", "r14", directed=False)
pgm.add_edge("r23", "r24", directed=False)
pgm.add_edge("r23", "r34", directed=False)
pgm.add_edge("r13", "r23", directed=False)
pgm.add_edge("r13", "r34", directed=False)
pgm.add_edge("r24", "r34", directed=False)
pgm.add_edge("r14", "r34", directed=False)
pgm.add_edge("r14", "r24", directed=False)

# Output
pgm.render()
pgm.figure.savefig("icml_2021_Fig.eps")