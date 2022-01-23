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


    void Test(Camera camera, Vector3 in_wp)
    {

        Vector4 wp = new Vector4(in_wp.x, in_wp.y, in_wp.z, 1.0f);
        Vector4 vp = camera.worldToCameraMatrix * wp;
        Vector4 pp = camera.projectionMatrix * vp;
        Vector3 hp = new Vector3(0.0f, 0.0f, 0.0f);
        if (pp.w != 0f) {
            hp.x = pp.x/pp.w;
            hp.y = pp.y/pp.w;
            hp.z = pp.w;
        }
        Vector3 sp = new Vector3((hp.x + 1f) * .5f * camera.pixelWidth, (hp.y + 1f) * .5f * camera.pixelHeight, hp.z);


        Vector3 wp_check = new Vector3(wp.x, wp.y, wp.z);
        Vector3 vp_check = camera.WorldToViewportPoint(wp_check);
        Vector3 sp_check = camera.WorldToScreenPoint(wp_check);
    }

    void EdgeTest(Camera camera, Vector3 wp0, float u0, Vector3 wp1, float u1) {
        Vector3 sp0 = camera.WorldToScreenPoint(wp0);
        Vector3 sp1 = camera.WorldToScreenPoint(wp1);
        
        float iw0 = 1.0f / sp0.z;
        float iw1 = 1.0f / sp1.z;
        float iu0 = iw0 * u0;
        float iu1 = iw1 * u1;

        float dx = sp1.x - sp0.x;
        float dy = sp1.y - sp0.y;
        float l = (float)System.Math.Sqrt((dx * dx) + (dy * dy));

        float diu = iu1 - iu0;
        float r1 = diu / l;
        float diu_dx = r1 * dx;
        float diu_dy = r1 * dy;

        float check = diu_dx * dx + diu_dy * dy;

        float diw = iw1 - iw0;
        float r2 = diw / l;
        float diw_dx = r2 * dx;
        float diw_dy = r2 * dy;
    }

    // Compute barycentric coordinates (u, v, w) for
    // point p with respect to triangle (a, b, c)
    Vector3 Barycentric(Vector2 p, Vector2 a,Vector2 b, Vector2 c)
    {
        Vector2 v0 = b - a, v1 = c - a, v2 = p - a;
        float den = v0.x * v1.y - v1.x * v0.y;
        float v = (v2.x * v1.y - v1.x * v2.y) / den;
        float w = (v0.x * v2.y - v2.x * v0.y) / den;
        float u = 1.0f - v - w;
        return new Vector3(u, v, w);
    }

    // Compute barycentric coordinates (u, v, w) for
    // point p with respect to triangle (a, b, c)
    Vector3 BarycentricV2(Vector2 p, Vector2 a, Vector2 b, Vector2 c)
    {
        Vector2 v0 = b - a, v1 = c - a, v2 = p - a;
        float d00 = Vector2.Dot(v0, v0);
        float d01 = Vector2.Dot(v0, v1);
        float d11 = Vector2.Dot(v1, v1);
        float d20 = Vector2.Dot(v2, v0);
        float d21 = Vector2.Dot(v2, v1);
        float denom = d00 * d11 - d01 * d01;
        float v = (d11 * d20 - d01 * d21) / denom;
        float w = (d00 * d21 - d01 * d20) / denom;
        float u = 1.0f - v - w;
        return new Vector3(u, v, w);
    } 

    Vector3 BarycentricV3(Vector2 p, Vector2 v1, Vector2 v2, Vector2 v3)
    {
        float w1 = ((v2.y - v3.y) * (p.x - v3.x) + (v3.x - v2.x) * (p.y - v3.y)) /
                ((v2.y - v3.y) * (v1.x - v3.x) + (v3.x - v2.x) * (v1.y - v3.y));
        float w2 = ((v3.y - v1.y) * (p.x - v3.x) + (v1.x - v3.x) * (p.y - v3.y)) /
                ((v2.y - v3.y) * (v1.x - v3.x) + (v3.x - v2.x) * (v1.y - v3.y));
        float w3 = 1.0f - w1 - w2;
        return new Vector3(w1, w2, w3);
    }

    void TestBarycentric() {
        Camera camera = GetComponent<Camera>();

        Vector3 v0 = new Vector3(0.0f, 0.0f, 0.0f);
        Vector3 v1 = new Vector3(5.0f, 0.0f, 0.0f);
        Vector3 v2 = new Vector3(5.0f, 5.0f, 0.0f);
        Vector3 sp0 = camera.WorldToScreenPoint(v0); 
        Vector3 sp1 = camera.WorldToScreenPoint(v1); 
        Vector3 sp2 = camera.WorldToScreenPoint(v2); 

        Vector2 sp0_v2 = new Vector2(sp0.x, sp0.y);
        Vector2 sp1_v2 = new Vector2(sp1.x, sp1.y);
        Vector2 sp2_v2 = new Vector2(sp2.x, sp2.y);


        Vector3 p = new Vector3(2.5f, 2.5f, 0.0f);
        Vector3 spp = camera.WorldToScreenPoint(p);
        Vector2 spp_v2 = new Vector2(spp.x, spp.y);
        Vector3 w0 = Barycentric(spp_v2, sp0_v2, sp1_v2, sp2_v2);
        Vector3 w1 = BarycentricV2(spp_v2, sp0_v2, sp1_v2, sp2_v2);
        Vector3 w2 = BarycentricV3(spp_v2, sp0_v2, sp1_v2, sp2_v2);

    }


    // Update is called once per frame
    void Test0()
    {
        Camera camera = GetComponent<Camera>();

        Vector3 v0 = new Vector3(0.0f, 0.0f, 0.0f);
        Vector3 v1 = new Vector3(5.0f, 0.0f, 0.0f);
        Vector3 v2 = new Vector3(5.0f, 5.0f, 0.0f);
        Test(camera, v0);
        Test(camera, v1);
        Test(camera, v2);
        Vector3 sp0 = camera.WorldToScreenPoint(v0); 
        Vector3 sp1 = camera.WorldToScreenPoint(v1); 
        Vector3 sp2 = camera.WorldToScreenPoint(v2); 

        EdgeTest(camera, v0, 0.0f, v2, 1.0f);

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
                                Vector3 vx0 = vertices[i0];
                                Vector3 vx1 = vertices[i1];
                                Vector3 wv0 = localToWorldMatrix.MultiplyPoint(vx0);
                                Vector3 wv1 = localToWorldMatrix.MultiplyPoint(vx1);
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

    void Update()
    {
        TestBarycentric ();
    }
}
