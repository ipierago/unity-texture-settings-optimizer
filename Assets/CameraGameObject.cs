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


    void CHECK(bool exp) {
        
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
    Vector3 BarycentricV1(Vector2 p, Vector2 a,Vector2 b, Vector2 c)
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
        Vector3 w0 = BarycentricV1(spp_v2, sp0_v2, sp1_v2, sp2_v2);
        Vector3 w1 = BarycentricV2(spp_v2, sp0_v2, sp1_v2, sp2_v2);
        Vector3 w2 = BarycentricV3(spp_v2, sp0_v2, sp1_v2, sp2_v2);
    }

    Vector3 ComputeBarycentricCoordsV4(float px, float py, float v1x, float v1y, float v2x, float v2y, float v3x, float v3y)
    {
        float w1 = ((v2y - v3y) * (px - v3x) + (v3x - v2x) * (py - v3y)) /
                ((v2y - v3y) * (v1x - v3x) + (v3x - v2x) * (v1y - v3y));
        float w2 = ((v3y - v1y) * (px - v3x) + (v1x - v3x) * (py - v3y)) /
                ((v2y - v3y) * (v1x - v3x) + (v3x - v2x) * (v1y - v3y));
        float w3 = 1.0f - w1 - w2;
        return new Vector3(w1, w2, w3);
    }

    float[] ComputeBarycentricCoords(float px, float py, float v1x, float v1y, float v2x, float v2y, float v3x, float v3y)
    {
        float w1 = ((v2y - v3y) * (px - v3x) + (v3x - v2x) * (py - v3y)) /
                ((v2y - v3y) * (v1x - v3x) + (v3x - v2x) * (v1y - v3y));
        float w2 = ((v3y - v1y) * (px - v3x) + (v1x - v3x) * (py - v3y)) /
                ((v2y - v3y) * (v1x - v3x) + (v3x - v2x) * (v1y - v3y));
        float w3 = 1.0f - w1 - w2;
        return new float[] {w1, w2, w3};
    }    

    void TestInterpolation() {
        Camera camera = GetComponent<Camera>();

        Vector3 v0 = new Vector3(0.0f, 0.0f, 0.0f);
        Vector2 uv0 = new Vector3(0.0f, 0.0f);
        Vector3 v1 = new Vector3(5.0f, 0.0f, 100.0f);
        Vector2 uv1 = new Vector3(1.0f, 0.0f);
        Vector3 v2 = new Vector3(5.0f, 5.0f, 100.0f);
        Vector2 uv2 = new Vector3(1.0f, 1.0f);

        Vector3 sp0 = camera.WorldToScreenPoint(v0); 
        Vector3 sp1 = camera.WorldToScreenPoint(v1); 
        Vector3 sp2 = camera.WorldToScreenPoint(v2); 

        // interpolant for 1/w
        float iw0 = 1.0f / sp0.z;
        float iw1 = 1.0f / sp1.z;
        float iw2 = 1.0f / sp2.z;

        // interpolant for u/w and v/w
        Vector2 iuv0 = new Vector2(uv0.x * iw0, uv0.y * iw0);
        Vector2 iuv1 = new Vector2(uv1.x * iw1, uv1.y * iw1);
        Vector2 iuv2 = new Vector2(uv2.x * iw2, uv2.y * iw2);

        const int numSteps = 11;
        Vector2 delta = new Vector2((sp1.x - sp0.x) / (numSteps - 1), (sp1.y - sp0.y) / (numSteps - 1));
        Vector2[] spArray = new Vector2[numSteps];
        Vector3[] wArray = new Vector3[numSteps];
        float[] iwArray = new float[numSteps];
        Vector2[] iuvArray = new Vector2[numSteps];
        Vector2[] uvArray = new Vector2[numSteps];
        for (int i = 0; i < numSteps; ++i) {
            Vector2 sp = new Vector2(sp0.x + i * delta.x, sp0.y + i * delta.y);
            Vector3 w = ComputeBarycentricCoordsV4(sp.x, sp.y, sp0.x, sp0.y, sp1.x, sp1.y, sp2.x, sp2.y);
            float iw = w.x * iw0 + w.y * iw1 + w.z * iw2;
            Vector2 iuv = w.x * iuv0 + w.y * iuv1 + w.z * iuv2;
            Vector2 uv = iuv / iw;
            spArray[i] = sp;
            wArray[i] = w;
            iwArray[i] = iw;
            iuvArray[i] = iuv;
            uvArray[i] = uv;
        }        
    }

    Vector2 ComputeInterpolatedUV(float sx, float sy, Vector3[] sxsywArray, float[] rhwArray, Vector2[] rhwuvArray) {
        float[] w = ComputeBarycentricCoords(sx, sy, sxsywArray[0].x, sxsywArray[0].y, sxsywArray[1].x, sxsywArray[1].y, sxsywArray[2].x, sxsywArray[2].y);
        float rhw = w[0] * rhwArray[0] + w[1] * rhwArray[1] + w[2] * rhwArray[2];
        Vector2 rhwuv = w[0] * rhwuvArray[0] + w[1] * rhwuvArray[1] + w[2] * rhwuvArray[2];
        Vector2 uv = rhwuv / rhw;
        return uv;
    }

    float[] ComputeMipUVAreasForTriangle(Vector3[] sxsywArray, Vector2[] uvArray) {
        
        // Compute screen positions and interpolants at each vertex
        float[] rhwArray = new float[3]; // 1/w interpolant
        Vector2[] rhwuvArray = new Vector2[3]; // uv/ w interpolant
        for (int i = 0; i < 3; ++i) {
            rhwArray[i] = 1.0f / sxsywArray[i].z;
            rhwuvArray[i] = rhwArray[i] * uvArray[i];
        }

        // Output array
        float[] uvAreasArray = new float[3];

        const float d = 1.0f; // Distance from each point to interpolate
        for (int i = 0; i < 3; ++i) {
            // Data for this point
            float sx = sxsywArray[i].x;
            float sy = sxsywArray[i].y;
            Vector2 uv = uvArray[i];

            // Compute interpolated UV for d pixel in both screen directions
            Vector2 uv_dx = ComputeInterpolatedUV(sx + d, sy, sxsywArray, rhwArray, rhwuvArray);
            Vector2 uv_dy = ComputeInterpolatedUV(sx, sy + d, sxsywArray, rhwArray, rhwuvArray);

            // Complete computation of partial derivatives
            Vector2 duv_dx = uv_dx - uv;
            Vector2 duv_dy = uv_dy - uv;

            // Compute uv area of parallelogram mapped by a single pixel
            float a_uv = Mathf.Abs(duv_dx.x * duv_dy.y - duv_dy.x * duv_dx.y) / (d * d);
            uvAreasArray[i] = a_uv;
        }
        return uvAreasArray;
    }

    void ComputeMipUVAreasForTriangle_Test() {
        const float k_dx = 256.0f; const float k_dy = 256.0f;
        for (int iScale = 1; iScale < 8; ++iScale) {
            float dx = k_dx / iScale; float dy = k_dy / iScale;
            Vector3[] sxsywArray = new Vector3[] {new Vector3(0.0f, 0.0f, 1.0f), new Vector3(dx, 0.0f, 1.0f), new Vector3(dx, dy, 1.0f)};
            Vector2[] uvArray = new Vector2[] {new Vector2(0.0f, 0.0f), new Vector2(1.0f, 0.0f), new Vector2(1.0f, 1.0f)};
            float[] uvAreaArray = ComputeMipUVAreasForTriangle(sxsywArray, uvArray);
            float[] mipLevelArray = new float[3];
            for (int i = 0; i < 3; ++i) {mipLevelArray[i] = uvAreaArray[i] * k_dx * k_dy;}
        }
    }

    float[] ComputeMipLevelsForMesh(Camera camera, Vector2 textureSizeInPixels, Vector3[] meshVertexArray, Vector2[] meshUVArray, int[] meshIndexArray) {
        int numTriangles = meshIndexArray.Length / 3;
        float[] meshMipLevelArray = new float[numTriangles * 3];
        for (int iTriangle = 0; iTriangle < numTriangles; ++iTriangle) {
            int i0 = meshIndexArray[iTriangle * 3];
            int i1 = meshIndexArray[iTriangle * 3 + 1];
            int i2 = meshIndexArray[iTriangle * 3 + 2];
            Vector3[] sxsywArray = new Vector3[] {camera.WorldToScreenPoint(meshVertexArray[i0]), camera.WorldToScreenPoint(meshVertexArray[i1]), camera.WorldToScreenPoint(meshVertexArray[i2])};
            Vector2[] uvArray = new Vector2[] {meshUVArray[i0], meshUVArray[i1], meshUVArray[i2]};
            float[] triangleMipUVAreasArray = ComputeMipUVAreasForTriangle(sxsywArray, uvArray);
            // Convert uv areas to texel areas which are equivalent to mip level
            for (int i = 0; i < 3; ++i) {meshMipLevelArray[iTriangle * 3 + i] = triangleMipUVAreasArray[i] * textureSizeInPixels.x * textureSizeInPixels.y;}
        }
        return meshMipLevelArray;
    }

    void TestMipSelection() {
        Camera camera = GetComponent<Camera>();

        Vector3[] meshVertexArray = new Vector3[] { new Vector3(0.0f, 0.0f, 0.0f), new Vector3(5.0f, 0.0f, 100.0f), new Vector3(5.0f, 5.0f, 100.0f) };
        Vector2[] meshUVArray = new Vector2[] {new Vector3(0.0f, 0.0f), new Vector3(1.0f, 0.0f), new Vector3(1.0f, 1.0f) };
        int[] meshIndexArray = new int[] {0, 1, 2};
        Vector2 textureSizeInPixels = new Vector2(256.0f, 256.0f);
        float [] meshMipLevelArray = ComputeMipLevelsForMesh(camera, textureSizeInPixels, meshVertexArray, meshUVArray, meshIndexArray);
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
        //TestBarycentric ();
        //TestInterpolation();
        //TestMipSelection();
        ComputeMipUVAreasForTriangle_Test();
    }
}
