
# pair covering graph

require(ggplot2)
require(network)
require(ggnet)

l1 = load("../DATA/PAIR/mat.list.pair.rda")
l2 = load("../OBJ/cover.pair.rda")
l3 = load("../DATA/PAIR/links.rda")

get_all_links = function(a){
  
  Reduce(rbind, lapply(a, 
                       function(b) t(sapply(1:(length(b)-1), 
                                            function(j) c(b[j], 
                                                          b[j+1], 
                                                          ifelse(j==1, -1, ifelse(j==length(b)-1, 1, 0))
                                            )
                       ) ))
  )
  
}

# we select breast as an example here

mat.sol = list.cover.pair$Breast$mat.sol
x = list.cover.pair$Breast$sol

u = links[x]

v = get_all_links(u)
x.s = unique(unname(sapply(x, function(a) unlist(strsplit(a, split="_"))[1] )))
x.t = unique(unname(sapply(x, function(a) unlist(strsplit(a, split="_"))[2] )))

v.s = unique(v[which(v[, 3] == "-1"), 1])
v.t = unique(v[which(v[, 3] == "1"), 2])

nn = network(v[, 1:2])
vertices = network.vertex.names(nn)
edge.sizes = rep(1, nrow(v))
edge.sizes[which(v[, 3] == "1")] = 2

cols = sapply(vertices, function(w){
  if(w %in% x.s && w %in% x.t)
    "maroon"
  else if(w %in% x.s)
    "orange" 
  else if(w %in% x.t)
    "cadetblue3" 
  else
    "darkolivegreen3"
})

labs = sapply(vertices, function(w){
  if(w %in% v.s && w %in% v.t)
    paste(w, "(S|T)", sep="")
  else if(w %in% v.s)
    paste(w, "(S)", sep="")
  else if(w %in% v.t)
    paste(w, "(T)", sep="")
  else
    w
})

lab.cols = sapply(vertices, function(w){
  if(w %in% x.s && w %in% x.t)
    "indianred4"
  else if(w %in% x.s)
    "darkorange3" #"violetred2" #
  else if(w %in% x.t)
    "steelblue" #"steelblue" #
  else
    "darkgreen"
})

alphas = rep(0.5, length(vertices))

png("breast_network.png", width=600, height=600)
tryCatch({
set.seed(1)
ggnet2(nn, label=TRUE, 
       node.color=cols,
       size=20,
       label.size = 6,
       label.color = "black", 
       edge.color = as.character(factor(edge.sizes, levels=c(1, 2), labels=c("gray70", "gray50"))),
       alpha=alphas,
       edge.size=edge.sizes)+
  theme(panel.background = element_rect(fill = "gray90", colour = "gray90"))
}, error = function(e){})
dev.off()






