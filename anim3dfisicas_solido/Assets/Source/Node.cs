using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Node
{
    public Vector3 position;
    public Vector3 velocity;
    public Vector3 acceleration;
    public Vector3 force;
    public float mass;
    public bool is_fixed;

    public Node()
    {
        this.Init();
    }

    public Node(Vector3 position, float mass = 1, bool is_fixed = false)
    {
        this.Init();
        this.position = position;
        this.mass = mass;
        this.is_fixed = is_fixed;
    }

    private void Init()
    {
        this.position = Vector3.zero;
        this.velocity = Vector3.zero;
        this.acceleration = Vector3.zero;
        this.force = Vector3.zero;
        this.mass = 0;
        this.is_fixed = false;
    }

    //compute forces with damping!
    public void ComputeForce(float damping)
    {
        //raw force would be F = m*a
        this.force += this.mass * this.acceleration;

        //but we also add damping, so we need to subtract the damping force from the current force.
        this.force -= this.mass * this.velocity * damping;
    }

    public void AddForceVector(Vector3 force)
    {
        this.force += force;
    }
}
