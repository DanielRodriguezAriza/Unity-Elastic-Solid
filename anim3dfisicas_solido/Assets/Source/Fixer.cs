using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Fixer : MonoBehaviour
{
    #region InEditorVariables
    public List<PhysicsManagerObject> fixable_objects;
    #endregion

    private struct FixedNode
    {
        public Node node;
        public Vector3 offset;
        public FixedNode(Node node, Vector3 offset)
        {
            this.node = node;
            this.offset = offset;
        }
    };

    private List<FixedNode> fixed_nodes;
    private Collider fixer_collider;
    private Vector3 old_position;

    // Start is called before the first frame update
    void Start()
    {   
        this.fixed_nodes = new List<FixedNode>();
        if (this.fixable_objects == null)
            return;

        fixer_collider = GetComponent<Collider>();
        this.old_position = this.transform.position;
        
        foreach (var fixable in fixable_objects)
        {
            Node[] nodes = fixable.GetNodes();
            //print("aqui llego 1");
            foreach (var node in nodes)
            {
                //print("aqui llego 2");
                if (this.IsInBounds(node))
                {
                    //print("aqui llego 3");
                    node.is_fixed = true;
                    Vector3 offset = node.position - this.transform.position;
                    FixedNode fixed_node = new FixedNode(node, offset);
                    this.fixed_nodes.Add(fixed_node);
                }
            }
        }
    }

    // Update is called once per frame
    void Update()
    {
        if (this.transform.position != this.old_position)
        {
            Vector3 displacement = this.transform.position - this.old_position;
            foreach (var fixed_node in this.fixed_nodes)
            {
                fixed_node.node.position += displacement;
            }
            this.old_position = this.transform.position;
        }
    }

    private bool IsInBounds(Node node)
    {
        return this.fixer_collider.bounds.Contains(node.position);
    }
}
