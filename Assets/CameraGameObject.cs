using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CameraGameObject : MonoBehaviour
{
    public MeshFilter[] meshFilterArray;
    public MeshRenderer[] meshRendererArray;

    // Start is called before the first frame update
    void Start()
    {

        
    }

    // Update is called once per frame
    void Update()
    {
        Camera camera = GetComponent<Camera>();

        List<MeshFilter> listMeshFilter = new List<MeshFilter>();
        List<MeshRenderer> listMeshRenderer = new List<MeshRenderer>();
        GameObject[] allObjects = Object.FindObjectsOfType<GameObject>() ;
        foreach(GameObject gameObject in allObjects) {
            if (gameObject.activeInHierarchy) {
                MeshFilter meshFilter = gameObject.GetComponent<MeshFilter>();
                MeshRenderer meshRenderer = gameObject.GetComponent<MeshRenderer>();
                if (meshFilter != null && meshRenderer != null) {
                    listMeshFilter.Add(meshFilter);
                    listMeshRenderer.Add(meshRenderer);
                }
            }
        }
        meshFilterArray = listMeshFilter.ToArray();
        meshRendererArray = listMeshRenderer.ToArray();

        for (int meshIndex = 0; meshIndex < meshFilterArray.Length; ++meshIndex) {
            MeshFilter meshFilter = meshFilterArray[meshIndex];
            MeshRenderer meshRenderer = meshRendererArray[meshIndex];
            Mesh mesh = meshFilter.mesh;

            Matrix4x4 localToWorldMatrix = meshRenderer.localToWorldMatrix;
            Vector3[] normals = mesh.normals;
            Vector3[] vertices = mesh.vertices;
            Vector2[] uv = mesh.uv;
            Material[] materials = meshRenderer.materials;
            if (mesh.subMeshCount <= materials.Length) {
                for (int subMeshIndex = 0; subMeshIndex < mesh.subMeshCount; ++subMeshIndex) {
                    Material material = materials[subMeshIndex];
                    if (material.mainTexture != null) {
                        int[] triangles = mesh.GetTriangles(subMeshIndex);
                        int numTriangles = triangles.Length / 3;
                        for (int triangleIndex = 0; triangleIndex < numTriangles; ++triangleIndex) {
                            for (int edgeIndex = 0; edgeIndex < 3; ++edgeIndex) {
                                int i0 = triangles[triangleIndex * 3 + edgeIndex];
                                int i1 = triangles[triangleIndex * 3 + ((edgeIndex + 1) % 3)];
                                Vector3 v0 = vertices[i0];
                                Vector3 v1 = vertices[i1];
                                Vector3 wv0 = localToWorldMatrix.MultiplyPoint(v0);
                                Vector3 wv1 = localToWorldMatrix.MultiplyPoint(v1);
                                Vector3 vv0 = camera.WorldToViewportPoint(wv0);
                                Vector3 vv1 = camera.WorldToViewportPoint(wv1);
                                Vector3 sv0 = camera.WorldToScreenPoint(wv0);
                                Vector3 sv1 = camera.WorldToScreenPoint(wv1);
                                Vector2 uv0 = uv[i0];
                                Vector2 uv1 = uv[i1];

                            }
                            
                        }
                    }
                }
            }
            

            

        }

    }
}
