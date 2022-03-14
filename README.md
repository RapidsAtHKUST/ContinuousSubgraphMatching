# Continuous Subgraph Matching
## Introduction

We propose to model CSM as incremental view maintenance (IVM) to capture the design space of existing algorithms, and implement six representative CSM algorithms, including IncIsoMatch, SJ-Tree, Graphflow, IEDyn, TurboFlux, and SymBi, in a common framework based on IVM. 

We further conduct extensive experiments to evaluate the overall performance of competing algorithms as well as study the effectiveness of individual techniques to pinpoint the key factors leading to the performance differences. We obtain the following new insights into the performance: (1) existing algorithms start the search from an edge in the query graph that maps to an updated data edge, potentially leading to many invalid partial results; (2) all matching orders are based on simple heuristics, which appear ineffective at times; (3) index updates dominate the query time on some queries; and (4) the algorithm with constant delay enumeration bears significant index update cost. 

Consequently, no algorithm dominate the others in all cases. Therefore, we give a few recommendations based on our experiment results. In particular, the SymBi index is useful for sparse queries or long running queries. The matching orders of IEDyn and TurboFlux work well on tree queries, those of Graphflow on dense queries or when both query and data graphs are sparse, and otherwise, we recommend SymBi's matching orders. 

For the details, please refer to our VLDB'2022 paper "An In-Depth Study of Continuous Subgraph Matching" by [Xibo Sun](https://github.com/xibosun), [Dr. Shixuan Sun](https://shixuansun.github.io/), [Prof. Qiong Luo](https://cse.hkust.edu.hk/~luo/) and [Prof. Bingsheng He](https://www.comp.nus.edu.sg/~hebs/). If you have any further question, please feel free to contact us.

Please cite our paper if you use our source code.

- "Xibo Sun, Shixuan Sun, Qiong Luo, and Bingsheng He. An In-Depth Study of Continuous Subgraph Matching. VLDB 2022."

## Compile

Our framework requires c++17 and GCC 7.x (or later). One can compile the code by executing the following commands. 

```shell
make
```

## Execute

After a successful compilation, the binary file is created under the `build/` directory. One can execute CSM using the following command.

```shell
build/csm -q <query-graph-path> -d <data-graph-path> -u <update-stream-path> -a <algorithm>
```

where `<algorithm>` is chosen from `sj-tree`, `graphflow`, `iedyn`, `turboflux`, and `symbi`.

### Commandline Parameters

Other commandline parameters supported by the framework are listed in the following table.

| Command Line Parameters | Description                                                     | Valid Value      | Default Value |
|-------------------------|-----------------------------------------------------------------|------------------|---------------|
| --max-results           | The max number of results to be found on each update operation. | 0-4294967295     | 4294967295    |
| --time-limit            | Time limit for the incremental matching phase (in seconds).     | 0-4294967295     | 3600          |
| --report-initial        | Perform initial matching or not.                                | on/off           | on            |
| --initial-time-limit    | Time limit for the initial matching phase (in seconds).         | 0-4294967295     | 4294967295    |
| --print-prep            | Print preprocessing results or not.                             | on/off           | on            |
| --print-enum            | Print matches results or not.                                   | on/off           | off           |
| --homo                  | Enable subgraph homomorphism or not.                            | on/off           | off           |

For example, if one requires the framework (1) to return after finding the first result on each update operation; and (2) to spend at most 1 hour (3600 seconds) on the incremental matching, then the command should be

```shell
build/csm -q <query-graph-path> -d <data-graph-path> -u <update-stream-path> -a <algorithm> --max-results 1 --time-limit 3600
```

## Input File Format
Both the input query graph and data graph are vertex- and edge-labeled. Each vertex is represented by a distinct unsigned integer (from 0 to 4294967295). There is at most one edge between two arbitrary vertices. 

### Query Graph

Each line in the query graph file represent a vertex or an edge.

1. A vertex is represented by `v <vertex-id> <vertex-label>`;
2. An edge is represented by `e <vertex-id-1> <vertex-id-2> <edge-label>`.

Both the two endpoint vertices of an edge should present before the edge. For example, 

```
v 0 0
v 1 0
e 0 1 0
v 2 1
e 0 2 1
e 2 1 2
```

### Initial Data Graph

The initial data graph file has the same format as the query graph file.

### Graph Update Stream

Graph update stream is a collection of insertions and deletions of a vertex or an edge.

1. A vertex insertion is represented by `v <vertex-id> <vertex-label>`;
2. A vertex deletion is represented by `-v <vertex-id> <vertex-label>`;
3. An edge insertion is represented by `e <vertex-id-1> <vertex-id-2> <edge-label>`;
4. An edge deletion is represented by `-e <vertex-id-1> <vertex-id-2> <edge-label>`;

The vertex or edge to be deleted should exist in the graph, and the label should be the same as that in the graph. If an edge is inserted to the data graph, both its endpoints should exist. For example,

```
v 3 1
e 2 3 2
-v 2 1
-e 0 1 0
```

## Datasets and Querysets

The graph datasets and their corresponding querysets used in our paper can be downloaded [here](https://hkustconnect-my.sharepoint.com/:f:/g/personal/xsunax_connect_ust_hk/Et-cxVY7l5FCoZoKeDyMzmQBaCBn8ffbPFFQfIFOqGIodA?e=4vT3OI).
