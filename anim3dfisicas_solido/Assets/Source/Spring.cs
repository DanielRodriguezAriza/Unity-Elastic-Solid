using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Spring
{
    public Node node_a, node_b;
    public float length0, length; //this.length0 is the starting length of the spring, and this.length is the current length of the spring.
    public float stiffness;

    public Spring(Node node_a, Node node_b, float stiffness)
    {
        this.node_a = node_a;
        this.node_b = node_b;
        this.stiffness = stiffness;
        this.UpdateLength();
        this.length0 = this.length;
    }

    public void UpdateLength()
    {
        //obtain length of the vector going from node_a to node_b
        //this will be the spring's length
        this.length = (this.node_a.position - this.node_b.position).magnitude;
    }

    public void ComputeForce(float damping)
    {
        Vector3 u = this.node_a.position - this.node_b.position;
        u.Normalize();
        Vector3 force = -this.stiffness * (this.length - this.length0) * u;
        this.node_a.force += force;
        this.node_b.force -= force;

        //damping logic:
        Vector3 vel_vector = this.node_a.velocity - this.node_b.velocity;
        Vector3 force_damping = -damping * this.stiffness * Vector3.Dot(u, vel_vector) * u;
        this.node_a.force += force_damping;
        this.node_b.force -= force_damping;
    }

}
