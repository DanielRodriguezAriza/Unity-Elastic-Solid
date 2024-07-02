using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/*
Note:
    In the original implementation, the springs and nodes had a visual representation, which was any mesh we chose, for example, spheres and cylinders.
    In the current implementation, to achieve the same as before, so as to prevent having to make a custom PhysicsManagerWhatever class for each and every single
    type of mass spring body we want to simulate, we can either:
        1) Make use of an INode interface and ISpring interface for the nodes and springs lists.
        or        
        2) Create a custom class similar to MassSpringCloth where we can define the meshes that we want to use and their positions to simulate the cloth
    for now, this code is going for option 2 in the case that I ever want to implement again a similar thing to what we originally had, but I have a hunch that this
    isn't the cleanest way to do things...
    We could also have a system where we can add some sort of "NodeComponent" and "SpringComponent" classes lists and iterate over them as well in the case that we allow
    (with a bool checkbox in the editor) some sort of "custom built spring mass shape" or whatever, so that we could retail the old system and iterate over
    those components and call the same functions, which, under the hood, will do the same calls as Node and Spring do.
    OR, instead of picking objects from the physics manager, we pick the physics manager from the objects... not sure about this one, but it could work.
*/

public class PhysicsManager : MonoBehaviour
{

    public enum Integration
    {
        EXPLICIT = 0,
        SYMPLECTIC = 1,
    };

    #region InEditorVariables
    public bool paused;
    public float time_step;
    public Vector3 gravity;
    public Integration integration_method;
    public List<Node> nodes;
    public List<Spring> springs;
    #endregion

    #region OtherVariables
    #endregion

    #region MonoBehaviour
    public PhysicsManager()
    { 
        
    }

    // Start is called before the first frame update
    void Start()
    {
        nodes = new List<Node>();
        springs = new List<Spring>();
    }

    // Update is called once per frame
    void Update()
    {
        if (Input.GetKeyUp(KeyCode.P))
            this.paused = !this.paused;
    }

    void FixedUpdate()
    {
        if (this.paused)
            return;

        switch (this.integration_method)
        {
            case Integration.EXPLICIT:
                this.stepExplicit();
                break;
            case Integration.SYMPLECTIC:
                this.stepSymplectic();
                break;
            default:
                throw new System.Exception("ERROR : THIS SHOULD NEVER HAPPEN, ELLIS, WTF DID YOU JUST DO????");
        }
    }

    #endregion

    #region StepFunctions
    private void stepExplicit()
    { }

    private void stepSymplectic()
    {
        foreach (var node in this.nodes)
        {
            node.force = Vector3.zero;
            node.ComputeForce(0);
        }

        foreach (var spring in this.springs)
        {
            spring.ComputeForce(0);
        }

        foreach (var node in this.nodes)
        {
            if (!node.is_fixed)
            {
                node.velocity += this.time_step / node.mass * node.force;
                node.position += time_step * node.velocity;
            }
        }

        foreach (var spring in this.springs)
        {
            spring.UpdateLength();
        }

    }

    #endregion

    #region OtherMethods
    public void AddNode(Node node)
    {
        this.nodes.Add(node);
    }

    public void AddSpring(Spring spring)
    {
        this.springs.Add(spring);
    }
    #endregion
}
