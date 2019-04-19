import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Scanner;
import java.util.Set;
import java.util.function.BiFunction;

/**
 * CSC 505 HW5 - Johnson's Algorithm on a Sparse Graph
 *
 * Reads a graph description and list of path queries from stdin, and returns
 * the lengths of the shortest paths for the queries.
 *
 */
public class johnsons {

    /**
     * Executes Johnson's Algorithm and prints the query path lengths.
     *
     * @param n number of vertices in the graph
     * @param m number of edges in the graph
     * @param edges set of edges of the graph
     * @param queries list of path queries
     */
    private void run(int n, int m, Set<Edge> edges, List<Query> queries) {
        Vertex[] vertices = new Vertex[n];
        for (Edge e : edges) {
            vertices[e.src.id] = e.src;
            vertices[e.dest.id] = e.dest;
        }

        Graph g = new Graph(vertices, edges);
        Integer[][] D = johnson(g);

        if (D != null) {
            for (Query query : queries) {
                Integer len = D[query.src][query.dest];
                System.out.println(String.format("%s -> %s = %s",
                        query.src,
                        query.dest,
                        len != null ? len : "x"));
            }
        }
    }

    /**
     * Johnson's Algorithm on a sparse graph. Returns an n x n matrix of
     * shortest path lengths between any two vertices.
     *
     * @param g Graph
     * @return an n x n matrix of shortest path lengths between any two vertices
     */
    private Integer[][] johnson(Graph g) {
        // new vertex not in G.V
        Vertex s = new Vertex(g.vertices.length, 0);

        // G.V  U  {s}
        Vertex[] gPrimeVertices = new Vertex[g.vertices.length + 1];
        for (Vertex v : g.vertices) {
            //gPrimeVertices[v.id] = new Vertex(v.id, v.distance);
            gPrimeVertices[v.id] = v;
        }
        gPrimeVertices[g.vertices.length] = s;

        // G.E  U  {(s, v) : v \in G.V} and w(s, v) = 0 for all v \in G.V
        Set<Edge> gPrimeEdges = new HashSet<>();
        for (Edge e : g.edges) {
            gPrimeEdges.add(new Edge(gPrimeVertices[e.src.id], gPrimeVertices[e.dest.id], e.weight));
        }

        for (Vertex v : gPrimeVertices) {
            gPrimeEdges.add(new Edge(s, v, 0));
        }

        Graph gPrime = new Graph(gPrimeVertices, gPrimeEdges);

        // Run Bellman-Ford to detect negative edge weight cycles
        if (!bellmanFord(gPrime, (u,v) -> weight(gPrime, u, v), s)) {
            System.out.println("Negative edge weight cycle");
            return null;
        }

        // For every vertex in G'.V, set h(v) = \delta(s, v) computed by Bellman-Ford
        Integer[] h = new Integer[gPrime.vertices.length];
        for (Vertex v : gPrime.vertices) {
            h[v.id] = v.distance;
        }

        Integer[][] D = new Integer[g.vertices.length][g.vertices.length];
        BiFunction<Vertex, Vertex, Integer> w = (u, v) -> weight(g, u, v);
        for (Vertex u : g.vertices) {
            // Run Dijkstra
            dijkstra(g, (x, y) -> reweight(w, x, y, h), u);
            for (Vertex v : g.vertices) {
                if (v.distance != Integer.MAX_VALUE) {
                    D[u.id][v.id] = v.distance + h[v.id] - h[u.id];
                }
            }
        }

        return D;
    }

    /**
     * Weight function to return the weight of an edge connecting vertices u and v.
     *
     * @param g Graph
     * @param u source vertex
     * @param v destination vertex
     * @return edge weight connecting vertices u and v.
     */
    private Integer weight(Graph g, Vertex u, Vertex v) {
        Set<Neighbor> neighbors = g.adjList.get(u.id);

        if (neighbors != null) {
            for (Neighbor neighbor : neighbors) {
                if (neighbor.vertexId == v.id) {
                    return neighbor.edgeWeight;
                }
            }
        }

        return null;
    }

    /**
     * Reweighting function to eliminate negative weights but preserve the
     * shortest length paths.
     *
     * @param w original edge weight function that could return negative weights
     * @param u source vertex
     * @param v destination vertex
     * @param h array of shortest path lengths from some starting vertex
     * @return new edge weight connecting vertices u, v
     */
    private Integer reweight(BiFunction<Vertex, Vertex, Integer> w, Vertex u, Vertex v, Integer[] h) {
        Integer weight = w.apply(u, v);
        if (weight == null) {
            return null;
        }

        return w.apply(u, v) + h[u.id] - h[v.id];
    }

    /**
     * Initializes the shortest path estimate and parent properties of every
     * vertex in the graph.
     *
     * @param g Graph
     * @param s start vertex.
     */
    private void initializeSingleSource(Graph g, Vertex s) {
        for (Vertex vertex : g.vertices) {
            vertex.distance = Integer.MAX_VALUE;
            vertex.parent = null;
        }

        s.distance = 0;
    }

    /**
     * Relaxes an edge connecting vertices u and v by calculating whether the
     * shortest from a start vertex to v can be improved by going through
     * vertex u.
     *
     * @param u source vertex
     * @param v destination vertex
     * @param w edge weight function
     */
    private void relax(Vertex u, Vertex v, BiFunction<Vertex, Vertex, Integer> w) {
        long relaxedDistance = ((long) u.distance) + w.apply(u, v);
        if (v.distance > relaxedDistance) {
            v.distance = (int) relaxedDistance;
            v.parent = u;
        }
    }

    /**
     * Bellman-Ford algorithm for shortest paths and negative weight cycle detection
     *
     * @param G Graph
     * @param w edge weight function
     * @param s start vertex
     * @return false if there is a negative weight cycle, true otherwise
     */
    private boolean bellmanFord(Graph G, BiFunction<Vertex, Vertex, Integer> w, Vertex s) {
        initializeSingleSource(G, s);

        for (int i = 1; i < G.vertices.length; i++) {
            for (Edge e : G.edges) {
                relax(e.src, e.dest, w);
            }
        }

        for (Edge e : G.edges) {
            if (e.dest.distance > ((long) e.src.distance) + w.apply(e.src, e.dest)) {
                return false;
            }
        }

        return true;
    }

    /**
     * Dijkstra's algorithm for shortest path length
     *
     * @param g Graph
     * @param w edge weight function
     * @param s start vertex
     */
    private void dijkstra(Graph g, BiFunction<Vertex, Vertex, Integer> w, Vertex s) {
        initializeSingleSource(g, s);

        List<Pair> pairs = new ArrayList<>();
        for(Vertex v : g.vertices) {
            pairs.add(new Pair(v.distance, v.id));
        }

        MinHeap queue = new MinHeap(2, g.vertices.length, pairs);

        while (queue.n != 0) {
            Vertex u = g.vertices[queue.removeMin().val];
            Set<Neighbor> neighbors = g.getAdjVertices(u.id);

            if (neighbors != null) {
                for (Neighbor neighbor : neighbors) {
                    relax(u, g.vertices[neighbor.vertexId], w);
                    queue.decreaseKey(neighbor.vertexId, g.vertices[neighbor.vertexId].distance);
                }
            }
        }
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        int n = scanner.hasNextInt() ? scanner.nextInt() : 0; // number of vertices
        int m = scanner.hasNextInt() ? scanner.nextInt() : 0; // number of edges

        Map<Integer, Vertex> vertices = new HashMap<>();
        List<Edge> edges = new LinkedList<>();
        for (int i = 0; i < m; i++) {
            Vertex src = vertices.computeIfAbsent(scanner.nextInt(), k -> new Vertex(k, Integer.MAX_VALUE));
            Vertex dest = vertices.computeIfAbsent(scanner.nextInt(), k -> new Vertex(k, Integer.MAX_VALUE));
            edges.add(new Edge(src, dest, scanner.nextInt()));
        }

        int k = scanner.nextInt();
        List<Query> queries = new LinkedList<>();
        for (int i = 0; i < k; i++) {
            queries.add(new Query(scanner.nextInt(), scanner.nextInt()));
        }

        new johnsons().run(n, m, new HashSet<>(edges), queries);
    }

    /**
     * Graph that uses an adjacency list representation
     *
     */
    static class Graph {
        private ArrayList<Set<Neighbor>> adjList; //Graph's adjacency list
        private Set<Edge> edges;
        private Vertex[] vertices; // Set of all the vertices in the graph

        /**
         * Initializes the Graph with given sets of vertices and edges
         *
         * @param vertices Given set of all vertices on the graph
         * @param edges set of all edges in the graph
         */
        private Graph(Vertex[] vertices, Set<Edge> edges) {
            this.vertices = vertices;
            this.edges = edges;
            this.adjList = new ArrayList<>(vertices.length);
            for (int i = 0; i < vertices.length; i++) {
                adjList.add(null);
            }

            for (Edge edge : edges) {
                Set<Neighbor> neighbors = adjList.size() > edge.src.id ? adjList.get(edge.src.id) : null;

                if (neighbors == null) {
                    neighbors = new HashSet<>();
                }

                neighbors.add(new Neighbor(edge.dest.id, edge.weight));
                adjList.set(edge.src.id, neighbors);
            }
        }

        /**
         * Returns adjacency list for a given vertex
         * @param u Vertex to get the adjacency list from
         * @return adjacency list of the given vertex
         */
        Set<Neighbor> getAdjVertices(int u) {
            return adjList.get(u);
        }

        /**
         * Returns an array of Vertex objects containing all the vertices on the graph
         * @return all the vertices of the graph
         */
        Vertex[] getVertices(){
            return this.vertices;
        }

        public Set<Edge> getEdges() {
            return edges;
        }
    }

    /**
     * Neighbor object
     * It represents a vertex on the adjacency list of another vertex from the graph
     */
    static class Neighbor {
        private int vertexId; //ID of this vertex
        private int edgeWeight; //Weight on the edge between this vertex and its parent

        /**
         * Initializes this object
         * @param vertexId ID for this vertex
         * @param edgeWeight Weight on the edge between this vertex and its parent
         */
        private Neighbor(int vertexId, int edgeWeight) {
            this.vertexId = vertexId;
            this.edgeWeight = edgeWeight;
        }

        /**
         * Neighbor hash code
         */
        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + edgeWeight;
            result = prime * result + vertexId;
            return result;
        }

        /**
         * Compares this neighbor with the given object
         * @param obj object to be compared
         * @return True if obj and this neighbor are the same
         */
        @Override
        public boolean equals(Object obj) {
            if (this == obj)
                return true;
            if (obj == null)
                return false;
            if (!(obj instanceof Neighbor))
                return false;
            Neighbor other = (Neighbor) obj;
            if (edgeWeight != other.edgeWeight)
                return false;
            if (vertexId != other.vertexId)
                return false;
            return true;
        }
    }

    static class Vertex {
        private int id;
        private int distance;
        private Vertex parent;

        public Vertex(int id, int distance) {
            this.id = id;
            this.distance = distance;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            Vertex vertex = (Vertex) o;
            return id == vertex.id;
        }

        @Override
        public int hashCode() {
            return Objects.hash(id);
        }
    }

    static class Edge {
        final Vertex src;
        final Vertex dest;
        final int weight;

        Edge(Vertex src, Vertex dest, int weight) {
            this.src = src;
            this.dest = dest;
            this.weight = weight;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            Edge edge = (Edge) o;
            return weight == edge.weight &&
                    src.equals(edge.src) &&
                    dest.equals(edge.dest);
        }

        @Override
        public int hashCode() {
            return Objects.hash(src, dest, weight);
        }
    }

    static class Query {
        final int src;
        final int dest;

        Query(int src, int dest) {
            this.src = src;
            this.dest = dest;
        }
    }

    /**
     * Key/value pair, used in the heap and its interface.
     * @author Dr. Sturgill
     */
    class Pair {
        public int key; //Key to organize the heap
        public int val; //Value of the heap node

        /**
         * Initializes a pair with key and value
         * @param k Given key
         * @param v Given value
         */
        public Pair( int k, int v ) {
            key = k;
            val = v;
        }
    }

    /**
     * Actual representation of the heap.
     *
     * @author Dr. Sturgill
     * @author Jeremy Maness
     */
    class MinHeap {
        private int p; // Power of 2 used as the branching factor
        private Pair[] tree; // Representation for the heap.
        private int n; // Number of elements in the heap.
        private int cap; // Capacity of the heap.
        private Integer[] vertexLocations; //array of vertices location to link vertices between the heap and the graph

        /**
         * Initializes the heap
         * @param p Power of 2 used as the branching factor
         * @param numVertices Maximum number of vertices that can be stored in the heap
         */
        public MinHeap( int p, int numVertices, List<Pair> pairs ) {
            this.p = p;
            cap = pairs.size(); //Initial capacity
            n = pairs.size(); //Number of vertices
            tree = pairs.toArray(new Pair[0]);
            vertexLocations = new Integer[numVertices];

            for (int i = 0; i < tree.length; i++) {
                vertexLocations[tree[i].val] = i;
            }

            buildMinHeap();
        }

        /**
         * Implements the Build-Min-Heap procedure from "Introduction to Algorithms, Third Edition"
         * (Cormen et al. 2009, p. 157).
         *
         * As shown on p. 157-159, a similar asymptotic analysis shows that this is O(V).
         *
         */
        private void buildMinHeap() {
            n = tree.length;
            for (int i = tree.length / 2; i >= 0; i--) {
                minHeapify(i);
            }
        }

        /**
         * Min heapify procedure that pushes the node at the specified index in
         * the heap down until the heap ordering constraint is satisfied.
         *
         * @param idx index of node in the heap
         */
        private void minHeapify(int idx) {

            /*
             * We need the branching factor below.
             */
            int branch = 1 << p;

            /*
             * Push this value down until it satisfies the ordering constraint.
             */
            int child = ( idx << p ) + 1;
            while ( child < n ) {
                // Find index of smallest child.
                int m = child;
                int end = child + branch;
                if ( end > n )
                    end = n;
                for ( int i = child + 1; i < end; i++ )
                    if (tree[i].key < tree[m].key)
                        m = i;

                /*
                 * Not happy about this early return.  Would be nice to have it in the condition
                 * on the loop.  Return early if we hit a point where we don't have to swap.
                 */
                if (tree[m].key >= tree[idx].key)
                    return;

                /*
                 * Swap the current value with its smallest child
                 */
                Pair temp = tree[ idx ];
                tree[ idx ] = tree[ m ];
                tree[ m ] = temp;

                vertexLocations[tree[idx].val] = idx;
                vertexLocations[tree[m].val] = m;

                /*
                 * Follow the value down into the tree.
                 */
                idx = m;
                child = ( idx << p ) + 1;
            }
        }

        /**
         * Remove the minimum value from the heap and executes a heapify operation to reorganize the heap.
         *
         * @return the vertex with minimum key value of the heap
         */
        Pair removeMin() {
            Pair v = tree[ 0 ];
            tree[ 0 ] = tree[ n - 1 ];
            n -= 1;

            vertexLocations[v.val] = null;
            vertexLocations[tree[0].val] = 0;

            minHeapify(0);

            return v;
        }

        /**
         * Decreases the key of the given vertex in the heap
         *
         * @param vertexId Given vertex ID
         * @param key Original vertex key value
         */
        void decreaseKey(int vertexId, int key) {
            Integer idx = vertexLocations[vertexId];

            if (idx == null) {
                return;
            }

            Pair pair = tree[idx];

            if (key > pair.key) {
                throw new RuntimeException("new key is larger than current key");
            }

            pair.key = key;

            int i = idx;
            while (i > 0 && tree[parent(i)].key > pair.key) {
                Pair temp = tree[i];
                tree[i] = tree[parent(i)];
                tree[parent(i)] = temp;

                vertexLocations[tree[parent(i)].val] = parent(i);
                vertexLocations[tree[i].val] = i;

                i = parent(i);
            }
        }

        /**
         * Returns the parent of the vertex with given index
         *
         * @param idx Vertex index
         * @return Index of the parent
         */
        int parent(int idx) {
            return (idx - 1) >> p;
        }
    }
}
