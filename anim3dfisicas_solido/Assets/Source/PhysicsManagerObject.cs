using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public abstract class PhysicsManagerObject : MonoBehaviour
{
    public abstract Node[] GetNodes();
    public abstract Spring[] GetSprings();
}
