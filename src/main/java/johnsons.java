import java.util.ArrayList;
import java.util.Arrays;
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
 *
 *
 */
public class johnsons {
    private int n;
    private int m;
    private Set<Edge> edges;
    private List<Query> queries;

    private Integer[][] shortestPaths;
    private Graph g;

    private johnsons(int n, int m, Set<Edge> edges, List<Query> queries) {
        this.n = n;
        this.m = m;
        this.edges = edges;
        this.queries = queries;
        this.shortestPaths = new Integer[n][n];

        Set<Vertex> vertices = new HashSet<>();
        for (Edge e : edges) {
            vertices.add(e.src);
            vertices.add(e.dest);
        }

        g = new Graph(vertices, edges);
    }

    private void run() {
        johnson(g);
    }

    private void johnson(Graph g) {
        Vertex s = new Vertex(g.vertices.length, 0);

        Set<Vertex> gPrimeVertices = new HashSet<>(Arrays.asList(g.vertices));
        gPrimeVertices.add(s);

        Set<Edge> gPrimeEdges = new HashSet<>(g.edges);
        for (Vertex v : g.vertices) {
            gPrimeEdges.add(new Edge(s, v, 0));
        }

        Graph gPrime = new Graph(gPrimeVertices, gPrimeEdges);

        // Run Bellman-Ford
        if (!bellmanFord(gPrime, (u,v) -> weight(gPrime, u, v), s)) {
            System.out.println("Negative edge weight cycle");
            return;
        }

        for (Vertex vertex : g.vertices) {
            // Run Dijkstra
            //dijkstra()
        }

        for (Query query : queries) {
            Integer len = getShortedPathLength(query.src, query.dest);

            System.out.println(String.format("%s -> %s = %s",
                    query.src,
                    query.dest,
                    len != null ? len : "x"));
        }
    }

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

    private boolean bellmanFord(Graph G, BiFunction<Vertex, Vertex, Integer> w, Vertex s) {
        initializeSingleSource(G, s);

        for (int i = 1; i < n; i++) {
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

    private void initializeSingleSource(Graph G, Vertex s) {
        for (Vertex vertex : G.vertices) {
            vertex.distance = Integer.MAX_VALUE;
            vertex.parent = null;
        }

        s.distance = 0;
    }

    private void relax(Vertex u, Vertex v, BiFunction<Vertex, Vertex, Integer> w) {
        long relaxedDistance = ((long) u.distance) + w.apply(u, v);
        if (v.distance > relaxedDistance) {
            v.distance = (int) relaxedDistance;
            v.parent = u;
        }
    }

    private Integer getShortedPathLength(int src, int dest) {
        return null; // TODO
    }


    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        int n = scanner.hasNextInt() ? scanner.nextInt() : 0;
        int m = scanner.hasNextInt() ? scanner.nextInt() : 0;

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

        new johnsons(n, m, new HashSet<>(edges), queries).run();
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
        private Graph(Set<Vertex> vertices, Set<Edge> edges) {
            this.vertices = new Vertex[vertices.size()];
            for (Vertex v : vertices) {
                this.vertices[v.id] = v;
            }

            this.edges = edges;
            this.adjList = new ArrayList<>(vertices.size());
            for (int i = 0; i < vertices.size(); i++) {
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
}
