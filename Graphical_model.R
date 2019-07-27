####グラフィカルモデルの作成####
library(DiagrammeR)

## LDAのグラフィカルモデル
grViz("
      digraph dot {
      graph [splines = line,compound = true, nodesep = .5, ranksep = .25, color = black, rankdir=LR]
      node [shape = circle,style = filled,fillcolor = white,color = black,label = '&alpha;'] alpha
      node [label = '&theta;@_{d}'] theta
      node [label = '&beta;'] beta
      node [label = '&phi;@_{k}'] phi
      node [label = 'z@_{dn}'] z
      node [fillcolor = grey,label = 'w@_{dn}'] w
      subgraph cluster1 {
      labelloc=b
      label = 'D = 1...d'
      edge [color = black]
      theta -> z
      subgraph cluster2 {
      labelloc=b
      label = 'n@_{d} = 1...N'
      edge [color = black]
      z -> w
      }
      }
      subgraph cluster3 {
      labelloc=b
      label = 'k = 1...K'
      edge [color = black]
      phi
      }
      
      edge [color = black]
      alpha -> theta
      beta -> phi -> w
      {rank = same; beta; w}
      }",engine = "dot")


## Hierarchical Structured Topic modelのグラフィカルモデル
grViz("
      digraph dot {
        graph [splines = line,compound = true, nodesep = .5, ranksep = .25, color = black, rankdir=LR]
        node [shape = circle,style = filled,fillcolor = white,color = black,label = '&alpha;'] alpha1
        node [label = '&alpha;@_{dt}'] alpha2
        node [label = '&theta;@_{d}'] theta1
        node [label = '&theta;@_{dt}'] theta2
        node [label = '&beta;@_{f}'] beta
        node [label = '&phi;@_{kf}@^{floor}'] phi
        node [label = 'z@_{dn}'] z
        node [fillcolor = grey,label = 'w@_{dn}'] w
        subgraph cluster1 {
          labelloc=b
          label = 'D = 1...d'
          edge [color = black]
          theta1 -> theta2
            subgraph cluster2 {
              labelloc=b
              label = 'day = 1...t'
              edge [color = black]
              alpha2 -> theta2
              {rank = same; alpha2; theta2}
          }
            subgraph cluster3 {
              labelloc=b
              label = 'n@_{d} = 1...N'
              edge [color = black]
              theta2 -> z -> w
              {rank = same; theta2; z; w}
          }
        }
        subgraph cluster4 {
          labelloc=b
          label = 'k = 1...K'
          edge [color = black]
          phi
        }
        
        edge [color = black]
        alpha1 -> theta1
        beta -> phi -> w
        {rank = same; beta; w}
      
      }",engine = "dot")


## Switching Hierarchical Structured Topic modelのグラフィカルモデル
grViz("
      digraph dot {
        graph [splines = line,compound = true, nodesep = .5, ranksep = .25, color = black, rankdir=TB, newrank=true]
        node [label = 'v2'] v2
        node [label = 'v1'] v1
        
        node [label = 'r@_{d}'] r
        node [label = '&alpha;@_{dt}'] alpha2
        node [label = '&theta;@_{d}'] theta1
        node [shape = circle,style = filled,fillcolor = white,color = black,label = '&alpha;'] alpha1
        node [label = '&theta;@_{dt}'] theta2
        node [label = '&gamma;'] gamma
        node [label = '&pi;@_{t}'] pi
        node [label = '&beta;@_{f,local}'] beta1
        node [label = '&beta;@_{f,global}'] beta2
        node [label = '&phi;@_{k@_{1}f}@^{local}'] phi1
        node [label = '&phi;@_{k@_{2}f}@^{global}'] phi2
        node [label = 'z@_{dn}@^{s@_{dn}}'] z
        node [label = 's@_{dn}'] s

        node [fillcolor = grey, label = 'w@_{dn}'] w
        subgraph cluster1 {
          labelloc=b
          label = 'D = 1...d'
          edge [color = black]
          theta1 -> theta2
          r -> s
            subgraph cluster2 {
              labelloc=b
              label = 'day = 1...t'
              edge [color = black]
              alpha2 -> theta2
              {rank = same; alpha2; theta2}
            }
            subgraph cluster3 {
              labelloc=b
              label = 'n@_{d} = 1...N'
              edge [color = black]
              theta2 -> z -> w
              s -> w
              {rank = same; theta2; z}
              {rank = same; s; w}
            }
        }
        subgraph cluster4 {
          labelloc=b
          label = 'k@_{2} = 1...K2'
          edge [color = black]
          phi2
        }
        subgraph cluster5 {
          labelloc=b
          label = 'k@_{1} = 1...K1'
          edge [color = black]
          phi1
        }
        subgraph cluster6 {
          labelloc=b
          label = 'day = 1...t'
          edge [color = black]
          pi
        }
        subgraph cluster7 {
          labelloc=b
          edge [color = black]
          v1 -> r
          v2 -> r
          
        }
        edge [color = black]
        gamma -> pi
        alpha1 -> theta1
        pi -> z
        beta1 -> phi1 -> w
        beta2 -> phi2 -> w
      {rank = same; v1; r}
      }",engine = "dot")
