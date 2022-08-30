using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CreatePrimitiveSample : MonoBehaviour
{
    public Material material;

    // Start is called before the first frame update
    void Start()
    {
        GameObject plane  = GameObject.CreatePrimitive(PrimitiveType.Plane);
        plane.GetComponent<Renderer>().material = material;
        plane.name = "MyPlane";

#if false
        GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);
        cube.transform.position = new Vector3(0, 0.5f, 0);
        cube.GetComponent<Renderer>().material = material;
#endif

#if false
        GameObject sphere = GameObject.CreatePrimitive(PrimitiveType.Sphere);
        sphere.transform.position = new Vector3(0, 1.5f, 0);
        sphere.GetComponent<Renderer>().material = material;
#endif

#if false
        GameObject capsule = GameObject.CreatePrimitive(PrimitiveType.Capsule);
        capsule.transform.position = new Vector3(2, 1, 0);
        capsule.GetComponent<Renderer>().material = material;
#endif

#if false
        GameObject cylinder = GameObject.CreatePrimitive(PrimitiveType.Cylinder);
        cylinder.transform.position = new Vector3(-2, 1, 0);
        cylinder.GetComponent<Renderer>().material = material;
#endif
    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
