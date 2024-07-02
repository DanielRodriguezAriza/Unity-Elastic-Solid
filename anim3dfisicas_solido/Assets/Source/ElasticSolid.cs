using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;
using UnityEngine.UIElements;
using static ElasticSolid;


public class ElasticSolid : PhysicsManagerObject
{

    #region Enums

    public enum StiffnessMode
    {
        StiffnessPerSpring = 0,
        StiffnessDensity
    }

    public enum MassMode
    {
        MassPerNode = 0,
        MassDensity
    }

    public enum Integration
    {
        EXPLICIT = 0,
        SYMPLECTIC = 1,
    };

    #endregion

    #region Structs

    // Stores edges as a struct where vertex_a is always smaller than vertex_b. This is done so that we can remove duplicates later on.
    public struct Edge
    {
        public int vertex_a;
        public int vertex_b;
        public Edge(int a, int b)
        {
            if (a <= b)
            {
                this.vertex_a = a;
                this.vertex_b = b;
            }
            else
            {
                this.vertex_a = b;
                this.vertex_b = a;
            }
        }
    };

    // Each Triangle defines a plane, defined by 3 points, that is, the triangle that composes each of the faces of the tetrahedra.
    public struct Triangle
    {
        public int[] nodes;

        public Triangle(int node1, int node2, int node3)
        {
            this.nodes = new int[3];
            this.nodes[0] = node1;
            this.nodes[1] = node2;
            this.nodes[2] = node3;
        }

        public void CalculateNormal(Node[] nodes, out Vector3 normal, out Vector3 faceCenterPos)
        {
            Vector3 node1pos = nodes[this.nodes[0]].position;
            Vector3 node2pos = nodes[this.nodes[1]].position;
            Vector3 node3pos = nodes[this.nodes[2]].position;
            faceCenterPos = (node1pos + node2pos + node3pos) / 3;
            normal = Vector3.Cross((node2pos - node1pos), (node3pos - node1pos)).normalized;
        }
    };

    // Each Tetrahedron defines a structure that contains 4 vertices and 4 planes (Triangles).
    // The nodes and faces are composed by maintaining the same ordering as defined by tetgen
    public struct Tetrahedron
    {
        #region Variables
        public int[] nodes;
        public Triangle[] triangles;
        public float volume;
        #endregion

        #region Constructors
        public Tetrahedron(int node1, int node2, int node3, int node4, Node[] nodes)
        {
            this.nodes = new int[4];
            this.triangles = new Triangle[4];
            this.volume = 0.0f;
            this.SetNodes(node1, node2, node3, node4);
            this.SetVolume(nodes);
        }
        #endregion

        #region PrivateMethods

        private void SetNodes(int node1, int node2, int node3, int node4)
        {
            this.nodes[0] = node1;
            this.nodes[1] = node2;
            this.nodes[2] = node3;
            this.nodes[3] = node4;
            this.UpdateTriangles();
        }

        private void UpdateTriangles()
        {
            this.triangles[0] = new Triangle(this.GetNode1(), this.GetNode2(), this.GetNode4()); // 1 2 4
            this.triangles[1] = new Triangle(this.GetNode4(), this.GetNode3(), this.GetNode1()); // 4 3 1
            this.triangles[2] = new Triangle(this.GetNode2(), this.GetNode3(), this.GetNode4()); // 2 3 4
            this.triangles[3] = new Triangle(this.GetNode1(), this.GetNode3(), this.GetNode2()); // 1 3 2
        }

        private void SetVolume(Node[] nodes)
        {
            Vector3 p0 = nodes[GetNode1()].position;
            Vector3 p1 = nodes[GetNode2()].position;
            Vector3 p2 = nodes[GetNode3()].position;
            Vector3 p3 = nodes[GetNode4()].position;
            this.volume = CalculateVolume(p0, p1, p2, p3, nodes);
        }

        #endregion

        #region PublicMethods

        public int GetNode1()
        {
            return this.nodes[0];
        }

        public int GetNode2()
        {
            return this.nodes[1];
        }

        public int GetNode3()
        {
            return this.nodes[2];
        }

        public int GetNode4()
        {
            return this.nodes[3];
        }

        public static float CalculateVolume(Vector3 p0, Vector3 p1, Vector3 p2, Vector3 p3, Node[] nodes)
        {
            return Mathf.Abs(Vector3.Dot((p1 - p0), Vector3.Cross((p2 - p0), (p3 - p0))));
        }

        public bool ContainsPoint(Vector3 point, Node[] nodes)
        {
            Vector3 n1, n2, n3, n4; // surface normals for each triangle of the tet
            Vector3 c1, c2, c3, c4; // center points for each triangle of the tet

            this.triangles[0].CalculateNormal(nodes, out n1, out c1);
            this.triangles[1].CalculateNormal(nodes, out n2, out c2);
            this.triangles[2].CalculateNormal(nodes, out n3, out c3);
            this.triangles[3].CalculateNormal(nodes, out n4, out c4);

            //signed distances
            float d1 = Vector3.Dot(n1, point - c1);
            float d2 = Vector3.Dot(n2, point - c2);
            float d3 = Vector3.Dot(n3, point - c3);
            float d4 = Vector3.Dot(n4, point - c4);

            // if all signs are the same, then the point is on the same side for all faces of the tet, which means that it is inside of it.
            bool sameSide = (d1 >= 0 && d2 >= 0 && d3 >= 0 && d4 >= 0) || (d1 <= 0 && d2 <= 0 && d3 <= 0 && d4 < 0);
            return sameSide;
        }

        #endregion
    }

    #endregion

    #region Classes

    // Class in charge of containing the data for the visual mesh that we're going to deform with the tet mesh.
    class MeshData
    {
        public Mesh mesh; //contains a ptr to the parent's mesh
        public Vector3[] vertices; //contains a list of the mesh's vertices
        public int[] triangles; //contains a list of the mesh's triangles
        public float[] weights; //contains a list of the weights for each vertex of the mesh relative to the tetahedra that contains it.
        public int[] containerTetrahedraIndex; // contains a list of the containing tets' indices.

        public MeshData(Mesh newMesh = null, Tetrahedron[] tetrahedra = null)
        {
            this.mesh = null;
            this.vertices = null;
            this.triangles = null;
            this.weights = null;
            this.containerTetrahedraIndex = null;
            CreateMeshData(newMesh, tetrahedra);
        }

        public void CreateMeshData(Mesh newMesh, Tetrahedron[] tetrahedra)
        {
            // Early return if any of the input params is null.
            if (newMesh == null || tetrahedra == null)
                return;

            // Set the mesh-related data
            this.mesh = newMesh; // GetComponentInChildren<MeshFilter>().mesh
            this.vertices = this.mesh.vertices;
            this.triangles = this.mesh.triangles;

            // Set the tet-related data.
            this.weights = new float[4 * vertices.Length];
            this.containerTetrahedraIndex = new int[this.vertices.Length];
        }

    }

    #endregion

    private Vector3 GetNormalVectorFromNodes(int node1, int node2, int node3)
    {
        Vector3 normal = Vector3.Cross((nodes[node2].position - nodes[node1].position), (nodes[node3].position - nodes[node1].position)).normalized;
        return normal;
    }

    private Vector3 GetCenterPointFromNodes(int node1, int node2, int node3)
    {
        Vector3 face_center_pos = (nodes[node1].position + nodes[node2].position + nodes[node3].position) / 3;
        return face_center_pos;
    }



    #region Variables
    //File variables (files generated with tetgen)
    [Header("TetGen Files")]
    [SerializeField] private TextAsset fileNode; // file that contains the nodes.
    [SerializeField] private TextAsset fileEle;  // file that contains the tetrahedra

    //Physics management public variables:
    [Header("Physics Simulation")]
    [SerializeField] private bool paused;
    [SerializeField] private float timeStep;
    [SerializeField] private Vector3 gravity;
    [SerializeField] private Integration integrationMethod; // determines the physics integration method used.
    [SerializeField] private int substeps; //determines the number of substeps to be performed for each FixedStep call.

    //Spring Cloth public variables:
    [Header("Node Config")]
    [SerializeField] private MassMode massMode = MassMode.MassPerNode;
    [SerializeField] private float nodeMass; // Contains a default value for the mass of each node.
    [SerializeField] private float nodeMassDensity;

    [Header("Spring Config")]
    [SerializeField] private StiffnessMode stiffnessMode = StiffnessMode.StiffnessPerSpring;
    [SerializeField] private float springStiffness; // Contains a default value for the stiffness of each spring.
    [SerializeField] private float springStiffnessDensity;

    [Header("Collider Config")]
    [SerializeField] private SphereCollider[] collidersPenaltyForces;
    [SerializeField] private SphereCollider[] collidersPullForces; //contains a list of the collider objects.

    //damping variables:
    [Header("Damping Config")]
    [SerializeField] private float dampingNode = 0.2f;
    [SerializeField] private float dampingSpring = 0.1f;

    //wind variables
    [Header("Wind Config")]
    [SerializeField] private Vector3 windVector; //this vector should always be normalized, so we'll normalize it on the Init() function.
    [SerializeField] private float windForce;

    //penalty force variables
    [Header("Penalty Forces Config")]
    [SerializeField] private float penaltyForce; //force applied due to penalty when colliding with sphere colliders.
    [SerializeField] private float penaltyPenetrationDistanceMargin; //margin of error of radius that must be added to smooth out the calculations.

    //pull force variables
    [Header("Pull Forces Config")]
    [SerializeField] private float pullForce;
    [SerializeField] private float pullPenetrationDistanceMargin;

    //other public variables:
    [Header("Debug Config")]
    [SerializeField] private bool debugEnabled; // determine if debug information is enabled or not
    [SerializeField] private bool debugPrint; // determines if debug information should be printed or not
    [SerializeField] private bool debugDraw; //determines if it should perform debug drawing during runtime.
    [SerializeField] private bool debugDrawObjectCenter; // draws a sphere in the point where the center of the object is located.
    [SerializeField] private bool debugDrawNodes; // draws a sphere for every node.
    [SerializeField] private bool debugDrawSprings; // draws a line for every spring.
    [SerializeField] private bool debugDrawNormals; // draws a dot in the center of each face and an outgoing line that represents the normal of the surface.
    #endregion

    #region OtherVariables
    //Nodes and Springs:
    private Node[] nodes;
    private Spring[] springs;
    private Tetrahedron[] tetrahedra;

    //Spring Cloth private variables:
    MeshData meshData;

    #endregion

    #region Constructors
    public ElasticSolid()
    { }
    #endregion

    #region Initialization
    public void Init()
    {
        // some empty initialization just in case the files do not load any data.
        this.nodes = new Node[0];
        this.springs = new Spring[0];
        this.tetrahedra = new Tetrahedron[0];

        this.meshData = new MeshData();

        this.windVector.Normalize(); //normalize wind vector just in case someone inputs a non-normalized wind vector!


        // First, load all of the data for the tetgen tetahedra.
        this.CreateNodes();
        this.CreateTetrahedra();
        this.CreateSprings();

        // Update Densities (mass density for nodes and stiffness density for springs):
        this.UpdateMassDensity();
        this.UpdateStiffnessDensity();

        // Then, generate the mesh data (requires the other step before so that it can make arrays of the proper length).
        this.CreateMeshData();

        // Calculate what the containing tets are.
        this.CreateContainingTetData();

        // Calculate weights for the mesh.
        this.CreateWeights();
    }
    #endregion

    #region MonoBehaviour
    // Start is called before the first frame update
    void Start()
    {

    }

    void Awake()
    {
        this.Init();
    }

    // Update is called once per frame
    void Update()
    {
        //print("Update");
        if (Input.GetKeyUp(KeyCode.P))
            this.paused = !this.paused;

        
        this.UpdateVertexPositions();
        this.UpdateSurfaceNormals();

    }

    void FixedUpdate()
    {
        if (this.paused)
            return;
        this.SubStep(this.substeps);
    }

    #endregion

    #region StepFunctions
    private void SubStep(int num)
    {
        for (int i = 0; i < num; ++i)
        {
            this.Step();
        }
    }
    private void Step()
    {
        switch (this.integrationMethod)
        {
            case Integration.EXPLICIT:
                this.stepExplicit();
                break;
            case Integration.SYMPLECTIC:
                this.stepSymplectic();
                break;
            default:
                throw new System.Exception("ERROR : THIS SHOULD NEVER HAPPEN, ELLIS, WTF DID YOU JUST DO???? EPISODE 2 : RETURN OF THE JEDI.");
        }
    }

    void SimulateForces()
    {
        this.SimulateWindForces();
        this.SimulatePenaltyForces();
        this.SimulatePullForces();
    }

    private void stepExplicit()
    {
        foreach (var node in this.nodes)
        {
            node.force = Vector3.zero;
            node.acceleration = this.gravity;
            node.ComputeForce(this.dampingNode);
        }
        this.SimulateForces();

        foreach (var spring in this.springs)
        {
            spring.ComputeForce(this.dampingSpring);
        }

        foreach (var node in this.nodes)
        {
            if (!node.is_fixed)
            {
                node.position += timeStep * node.velocity;
                node.velocity += this.timeStep / node.mass * node.force;
            }
        }

        foreach (var spring in this.springs)
        {
            spring.UpdateLength();
        }
    }

    private void stepSymplectic()
    {
        foreach (var node in this.nodes)
        {
            node.force = Vector3.zero;
            node.acceleration = this.gravity;
            node.ComputeForce(this.dampingNode);
        }
        this.SimulateForces();

        foreach (var spring in this.springs)
        {
            spring.ComputeForce(this.dampingSpring);
        }

        foreach (var node in this.nodes)
        {
            if (!node.is_fixed)
            {
                node.velocity += this.timeStep / node.mass * node.force;
                node.position += timeStep * node.velocity;
            }
        }

        foreach (var spring in this.springs)
        {
            spring.UpdateLength();
        }

    }

    #endregion

    #region OtherMethods
    #endregion

    #region ElasticSolidMethods
    private void CreateMeshData()
    {
        var msh = GetComponentInChildren<MeshFilter>().mesh;
        this.meshData.CreateMeshData(msh, tetrahedra);
    }

    private void CreateNodes()
    {
        if (this.fileNode != null)
            this.nodes = TGParser.ParseFileNode(this.fileNode, this.transform);
    }

    private void CreateSprings()
    {
        Edge[] edges = this.CreateEdgesFromTetrahedra();
        this.CreateSpringsFromEdges(edges);
    }

    // connect the 4 nodes of the tetahedron a,b,c,d
    // unique connections within this tet are: (6 edges in total)
    /*
        a,b : 1,2
        a,c : 1,3
        a,d : 1,4
        b,c : 2,3
        b,d : 2,4
        c,d : 3,4
    */
    private Edge[] CreateEdgesFromTetrahedra()
    {
        HashSet<Edge> edges = new HashSet<Edge>();
        foreach (var tetrahedron in this.tetrahedra)
        {
            Edge e1 = new Edge(tetrahedron.GetNode1(), tetrahedron.GetNode2());
            Edge e2 = new Edge(tetrahedron.GetNode1(), tetrahedron.GetNode3());
            Edge e3 = new Edge(tetrahedron.GetNode1(), tetrahedron.GetNode4());
            Edge e4 = new Edge(tetrahedron.GetNode2(), tetrahedron.GetNode3());
            Edge e5 = new Edge(tetrahedron.GetNode2(), tetrahedron.GetNode4());
            Edge e6 = new Edge(tetrahedron.GetNode3(), tetrahedron.GetNode4());
            edges.Add(e1);
            edges.Add(e2);
            edges.Add(e3);
            edges.Add(e4);
            edges.Add(e5);
            edges.Add(e6);
        }
        return edges.ToArray();
    }


    private void CreateSpringsFromEdges(Edge[] edges)
    {
        List<Spring> springs_list = new List<Spring>();
        foreach (var edge in edges)
        {
            int idx_node_a = edge.vertex_a;
            int idx_node_b = edge.vertex_b;
            Node node_a = this.nodes[idx_node_a];
            Node node_b = this.nodes[idx_node_b];
            Spring spring = new Spring(node_a, node_b, 0.0f); // Use a stiffness of 0 so that we can customize it through the addition loop later on when accounting for stiffness density.
            springs_list.Add(spring);
        }
        this.springs = springs_list.ToArray();
    }

    private void CreateTetrahedra()
    {
        if (this.fileEle != null)
            this.tetrahedra = TGParser.ParseFileEle(this.fileEle, this.transform, this.nodes);
    }

    private void CreateContainingTetData()
    {
        for (int i = 0; i < this.meshData.vertices.Length; ++i)
        {
            var vertex = this.meshData.vertices[i];
            vertex = this.transform.TransformPoint(vertex);
            for (int j = 0; j < tetrahedra.Length; ++j)
            {
                var tet = tetrahedra[j];
                if (tet.ContainsPoint(vertex, nodes))
                {
                    this.meshData.containerTetrahedraIndex[i] = j;
                    print($"Vertex {i} is contained by Tetrahedron {j}");
                    break;
                }
                else
                {
                    this.meshData.containerTetrahedraIndex[i] = -1; // indicates that no tet contains this vertex.
                }
            }
        }
    }

    private void CreateWeights()
    {
        for (int i = 0; i < meshData.vertices.Length; ++i)
        {
            int startIdx = i * 4;
            int tetIdx = meshData.containerTetrahedraIndex[i];
            Vector3 position = transform.TransformPoint(meshData.vertices[i]);

            // skip those vertices that are outside of the tet mesh
            if (tetIdx < 0)
            {
                continue;
            }

            if (tetIdx >= tetrahedra.Length)
            {
                continue;
            }

            Vector3 nodePos1 = nodes[tetrahedra[tetIdx].GetNode1()].position;
            Vector3 nodePos2 = nodes[tetrahedra[tetIdx].GetNode2()].position;
            Vector3 nodePos3 = nodes[tetrahedra[tetIdx].GetNode3()].position;
            Vector3 nodePos4 = nodes[tetrahedra[tetIdx].GetNode4()].position;

            meshData.weights[startIdx + 0] = Tetrahedron.CalculateVolume(nodePos2, nodePos3, position, nodePos4, nodes) / tetrahedra[tetIdx].volume;
            meshData.weights[startIdx + 1] = Tetrahedron.CalculateVolume(nodePos1, nodePos3, nodePos4, position, nodes) / tetrahedra[tetIdx].volume;
            meshData.weights[startIdx + 2] = Tetrahedron.CalculateVolume(nodePos1, nodePos2, position, nodePos4, nodes) / tetrahedra[tetIdx].volume;
            meshData.weights[startIdx + 3] = Tetrahedron.CalculateVolume(nodePos1, nodePos2, nodePos3, position, nodes) / tetrahedra[tetIdx].volume;
        }
    }


    //updates the vertex positions of the mesh based on the node's position data.
    private void UpdateVertexPositions()
    {
        for (int i = 0; i < this.meshData.vertices.Length; ++i)
        {
            // The position of each vertex is calculated as:
            // p = sum(wi * pi) for all of the 4 nodes of the tet that contains this point.
            // this.meshData.vertices[i] = this.transform.InverseTransformPoint(this.nodes[i].position);

            // If this vertex is not contained by any tet (caused by errors or by fixers), then we skip it:
            if (meshData.containerTetrahedraIndex[i] < 0)
            {
                continue;
            }

            float w1 = meshData.weights[i * 4 + 0];
            float w2 = meshData.weights[i * 4 + 1];
            float w3 = meshData.weights[i * 4 + 2];
            float w4 = meshData.weights[i * 4 + 3];

            var tet = tetrahedra[meshData.containerTetrahedraIndex[i]];
            Vector3 p1 = nodes[tet.GetNode1()].position;
            Vector3 p2 = nodes[tet.GetNode2()].position;
            Vector3 p3 = nodes[tet.GetNode3()].position;
            Vector3 p4 = nodes[tet.GetNode4()].position;

            Vector3 pos = w1 * p1 + w2 * p2 + w3 * p3 + w4 * p4;
            this.meshData.vertices[i] = this.transform.InverseTransformPoint(pos);
            //this.meshData.vertices[i] = this.meshData.vertices[i];
        }
        this.meshData.mesh.vertices = this.meshData.vertices;
    }

    private void UpdateSurfaceNormals()
    {
        if (meshData.mesh == null)
            return;
        
        this.meshData.mesh.RecalculateNormals();
        this.meshData.mesh.RecalculateTangents();
    }


    //Reminder: Wind force is calculated and applied over triangles rather than vertices so taht we can evenly distribute the wind force across the entire face surface to all vertices (for instance, 1/3 of the force applied to the surface is distributed to each vertex that makes up the triangle, etc...).
    //We could theoretically calculate it for every single vertex, but the reason we do so for every single triangle is because we think of the mesh as a surface, and the wind collides with multiple points at the same time, not with every single point individually, which means that the mesh is affected by the normal of the surface. Since we don't have infinitely many particles / points, we need to depend on surface normals for this. Or, we could store our own vertex normals on every single node, but we'll go with this implementation instead.
    //As it is currently implemented, the wind force is calculated to be a vector that goes in the direction of the surface normal. We calculate its force according to the wind force, and the angle between the wind vector and the surface normal (that's what the dot product is for!).
    private Vector3 CalculateWindForce(Node a, Node b, Node c, Vector3 wind_vector, float wind_force)
    {
        Vector3 v1 = b.position - a.position;
        Vector3 v2 = c.position - a.position;
        Vector3 cross_product = Vector3.Cross(v1, v2);

        Vector3 triangle_normal = cross_product.normalized;
        float triangle_area = cross_product.magnitude / 2;

        Vector3 ans = (Vector3.Dot(triangle_normal, wind_vector) * wind_force * triangle_area / 3) * triangle_normal;

        return ans;
    }

    //Extra notes about cross product:
    /*
        To calculate the normal of the triangle formed by nodes {A,B,C}, we can do the following: 
        
        Since we know that each surface is a flat triangle, we know that we can use 2 vectors that go from one point to another point each from the triangle to
        calculate the cross product and obtain a vector which will ne perpendicular to both of them, which we can normalize to obtain a valid normal vector.
        
        Steps to follow:
            1) calculate vector BA
            2) calculate vector CA
            3) cross product of BA and CA
            4) normalize the obtained vector and return it
        
        EZ.

        Other important info obtained from cross product is:
            -area of the polygon defined by the 2 vectors (if we divide this by 2, we get the area of a triangle. Again, ez).
    */

    private void SimulateWindForces()
    {
        foreach (var tet in this.tetrahedra)
        {
            foreach (var tri in tet.triangles)
            {
                var node0 = this.nodes[tri.nodes[0]];
                var node1 = this.nodes[tri.nodes[1]];
                var node2 = this.nodes[tri.nodes[2]];
                Vector3 wind_force_per_node = CalculateWindForce(node0, node1, node2, this.windVector, this.windForce);
                node0.AddForceVector(wind_force_per_node);
                node1.AddForceVector(wind_force_per_node);
                node2.AddForceVector(wind_force_per_node);
            }
        }
    }

    private void DebugDraw()
    {
        // early return if debug drawing is not enabled.
        if (!this.debugEnabled)
            return;

        // simple failsafe early return in case any of the members are not properly initialized.
        if (this.springs == null || this.nodes == null || this.tetrahedra == null)
            return;

        if (this.debugPrint)
            print($"Nodos: {this.nodes.Length}; Muelles: {this.springs.Length};");

        // early return if debug drawing is disabled completely since everything after this point of the function involves drawing stuff.
        if (!this.debugDraw)
            return;

        if (this.debugDrawObjectCenter)
        {
            Gizmos.color = Color.yellow;
            Gizmos.DrawSphere(transform.position, 0.2f);
        }

        if (this.debugDrawSprings)
        {
            Gizmos.color = Color.red;
            foreach (var spring in this.springs)
            {
                //Gizmos.color = Color.red; // we could put this outside of the loop for performance reasons I suppose...
                Gizmos.DrawLine(spring.node_a.position, spring.node_b.position);
            }
        }

        if (this.debugDrawNodes)
        {
            Gizmos.color = Color.blue;
            foreach (var node in this.nodes)
            {
                //Gizmos.color = Color.blue; // we could put this outside of the loop for performance reasons I suppose...
                Gizmos.DrawSphere(node.position, 0.2f);
            }
        }

        if (this.debugDrawNormals)
        {
            Gizmos.color = Color.cyan;
            foreach (var tet in this.tetrahedra)
            {
                foreach (var triangle in tet.triangles)
                {
                    Vector3 node1pos = this.nodes[triangle.nodes[0]].position;
                    Vector3 node2pos = this.nodes[triangle.nodes[1]].position;
                    Vector3 node3pos = this.nodes[triangle.nodes[2]].position;
                    Vector3 face_center_pos = (node1pos + node2pos + node3pos) / 3;
                    Vector3 normal = Vector3.Cross((node2pos - node1pos), (node3pos - node1pos)).normalized;
                    //Gizmos.color = Color.cyan;
                    Gizmos.DrawLine(face_center_pos, face_center_pos + normal);
                    //Gizmos.color = Color.cyan;
                    Gizmos.DrawSphere(face_center_pos, 0.04f);
                }
            }
        }
    }

    public void OnDrawGizmos()
    {
        this.DebugDraw();
    }

    private Vector3 CalculatePenetrationForce(SphereCollider collider, Vector3 point, float penetration_distance_margin, float penetration_force)
    {
        Vector3 force = Vector3.zero;

        Vector3 vector_director = collider.transform.position - point;

        float penetration_distance = vector_director.magnitude;

        float max_scale = Mathf.Max(collider.transform.localScale.x, collider.transform.localScale.y, collider.transform.localScale.z);

        if ((penetration_distance + penetration_distance_margin) < (collider.radius * max_scale))
        {
            Vector3 vector_normal = vector_director.normalized;
            force = penetration_force * penetration_distance * vector_normal;
        }
        return force;
    }

    private void CalculatePenaltyForce(SphereCollider collider, Node node)
    {
        Vector3 force = CalculatePenetrationForce(collider, node.position, this.penaltyPenetrationDistanceMargin, this.penaltyForce);
        node.force -= force;
    }

    private void CalculatePullForce(SphereCollider collider, Node node)
    {
        Vector3 force = CalculatePenetrationForce(collider, node.position, this.pullPenetrationDistanceMargin, this.pullForce);
        node.force += force;
    }

    private void SimulatePenaltyForces()
    {
        foreach (var node in this.nodes)
        {
            foreach (var collider in this.collidersPenaltyForces)
                this.CalculatePenaltyForce(collider, node);
        }
    }

    private void SimulatePullForces()
    {
        foreach (var node in this.nodes)
        {
            foreach (var collider in this.collidersPullForces)
                this.CalculatePullForce(collider, node);
        }
    }

    #endregion

    #region IFixableMethods
    public override Node[] GetNodes()
    {
        return this.nodes;
    }

    public override Spring[] GetSprings()
    {
        return this.springs;
    }

    #endregion

    #region PrivateMethods

    private void UpdateMassDensity()
    {
        switch(this.massMode)
        {
            default:
            case MassMode.MassPerNode:
                for (int i = 0; i < nodes.Length; ++i)
                    nodes[i].mass = nodeMass;
                break;

            case MassMode.MassDensity:
                foreach (var tetrahedron in tetrahedra)
                    foreach (var nodeIdx in tetrahedron.nodes)
                        nodes[nodeIdx].mass += nodeMassDensity * tetrahedron.volume * (1.0f / 4.0f); // 1/4 because there are 4 nodes in 1 tet.
                break;
        }
    }

    private void UpdateStiffnessDensity()
    {
        switch (this.stiffnessMode)
        {
            default:
            case StiffnessMode.StiffnessPerSpring:
                for (int i = 0; i < springs.Length; ++i)
                    springs[i].stiffness = springStiffness;
                break;

            case StiffnessMode.StiffnessDensity:
                for (int i = 0; i < springs.Length; ++i)
                    foreach (var tetrahedron in tetrahedra)
                        if ( // If the spring corresponds to this tetrahedron, then set the stiffness.
                            (
                            nodes[tetrahedron.GetNode1()] == springs[i].node_a ||
                            nodes[tetrahedron.GetNode2()] == springs[i].node_a ||
                            nodes[tetrahedron.GetNode3()] == springs[i].node_a ||
                            nodes[tetrahedron.GetNode4()] == springs[i].node_a
                            )
                            &&
                            (
                            nodes[tetrahedron.GetNode1()] == springs[i].node_b ||
                            nodes[tetrahedron.GetNode2()] == springs[i].node_b ||
                            nodes[tetrahedron.GetNode3()] == springs[i].node_b ||
                            nodes[tetrahedron.GetNode4()] == springs[i].node_b
                            )
                            )
                        {
                            springs[i].stiffness += (springStiffnessDensity * tetrahedron.volume / (springs[i].length0 * springs[i].length0)) * (1.0f / 6.0f); // 1/6 because there are 6 edges in 1 tet.
                            print($"stiffness is {springs[i].stiffness}");
                        }
                break;
        }
    }

    #endregion

}
