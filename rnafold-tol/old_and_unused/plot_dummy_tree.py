from ete3 import Tree, TreeStyle, faces, AttrFace, CircleFace, RectFace

# Loads a tree with internal node names
t = Tree("((((A:0.8,A:0.23),(A:0.3,A:0.12)),F):0.2,((D,((B:0.1,B:0.32),(B:0.24,B:0.15))):0.3,G:0.4):0.2);", format=1)

def layout(node):
    if node.is_leaf():
        if node.name=="A":
            node.color="Blue"
            node.attr2="Green"
        elif node.name=="B":
            node.color="Red"
            node.attr2="Orange"
        else:
            node.color="Grey"
            node.attr2="Yellow"
            
        faces.add_face_to_node( CircleFace(radius=5, color=node.color), node, column=1)
        faces.add_face_to_node( RectFace(width=10, height=10, fgcolor=node.attr2, bgcolor=node.attr2),   node, column=2)
        
ts = TreeStyle()
ts.layout_fn = layout
t.render("dummy_tree.png", tree_style=ts, h=700, w=800, dpi=500 )



