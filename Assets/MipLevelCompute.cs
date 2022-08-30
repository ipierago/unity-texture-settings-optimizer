﻿using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MipLevelCompute
{
   // Output for TransformStep
    struct TransformStepOutput {
        public Vector4[] hxhyhzwArray;
        public float[,] dotClipPlanesArray;
        public int [] clipCodeArray;
        public int clipCodeAND;
        public int clipCodeOR;
    }

    // Transform local positions to homogeneous space and compute dot products with clip planes and corresponding clip codes
    static TransformStepOutput TransformStep(Matrix4x4 worldViewProjMatrix, Vector4[] homoClipPlaneNormals, Vector3[] localPosArray) {
        TransformStepOutput output;

        // Allocate output arrays
        output.hxhyhzwArray = new Vector4[localPosArray.Length];
        output.dotClipPlanesArray = new float[localPosArray.Length, homoClipPlaneNormals.Length];
        output.clipCodeArray = new int[localPosArray.Length];
        output.clipCodeAND = 0;
        for (int i = 0; i < homoClipPlaneNormals.Length; ++i) output.clipCodeAND |= (1 << i);
        output.clipCodeOR = 0;

        for (int i = 0; i < localPosArray.Length; ++i) {
            // Homogeneous position
            Vector3 localPos3 = localPosArray[i];
            Vector4 localPos4 = new Vector4(localPos3.x, localPos3.y, localPos3.z, 1.0f);
            Vector4 homoPos = worldViewProjMatrix * localPos4;
            output.hxhyhzwArray[i] = homoPos;

            // Dot product with clip planes and corresponding clip code
            int clipCode = 0;
            for (int j = 0; j < homoClipPlaneNormals.Length; ++j) {
                float dot = Vector4.Dot(homoPos, homoClipPlaneNormals[j]);
                output.dotClipPlanesArray[i, j] = dot;
                if (dot < 0) clipCode |= (1 << j);
            }
            output.clipCodeArray[i] = clipCode;
            output.clipCodeAND &= clipCode;
            output.clipCodeOR |= clipCode;
        }
        return output;
    }

    struct ClipPolygonOutput {
        public List<Vector4> hxhyhzwList;
        public List<Vector2> uvList;
    }

    // Clip a polygon with UVs against the input homogeneous clip plane
    static ClipPolygonOutput ClipPolygon(List<Vector4> in_hxhyhzwList, List<Vector2> in_uvList, Vector4 in_clipPlaneNormal) {
        // Alloc output struct
        ClipPolygonOutput output;
        output.hxhyhzwList = new List<Vector4>();
        output.uvList = new List<Vector2>();

        // Walk each vertex and determine if the line from the previous vertex is clipped
        Vector4 prevPos = in_hxhyhzwList[in_hxhyhzwList.Count - 1];
        Vector2 prevUV = in_uvList[in_hxhyhzwList.Count - 1];
        float prevDot = Vector4.Dot(in_clipPlaneNormal, prevPos);
        for (int i = 0; i < in_hxhyhzwList.Count; ++i) {
            Vector4 thisPos = in_hxhyhzwList[i];
            Vector2 thisUV = in_uvList[i];
            float thisDot = Vector4.Dot(in_clipPlaneNormal, thisPos);
            if (prevDot * thisDot < 0) {
                // The edge is clipped.  Linearly interpolate to find the intersection point and add it to the new polygon
                float a = Mathf.Abs(prevDot) / (Mathf.Abs(prevDot) + Mathf.Abs(thisDot));
                Vector4 newPos = Vector4.Lerp(prevPos, thisPos, a);
                Vector2 newUV = Vector2.Lerp(prevUV, thisUV, a);
                output.hxhyhzwList.Add(newPos);
                output.uvList.Add(newUV);
            }
            // If this current vertex is not clipped, add it to our new polygon.
            if (thisDot > 0) {
                output.hxhyhzwList.Add(thisPos);
                output.uvList.Add(thisUV);
            }
            // Update prev values
            prevPos = thisPos;
            prevUV = thisUV;
            prevDot = thisDot;
        }

        return output;
    }

    // Output of clip step
    struct ClipStepOutput {
        public Vector4[] hxhyhzwArray;
        public Vector2[] uvArray;
        public int[] triListIndexArray;
    }

    // Clip input homogeneous geometry and UVs and return clipped geometry and UVs
    static ClipStepOutput ClipStep(Vector4[] in_hxhyhzwArray, Vector2[] in_uvArray, int[] in_triListIndexArray, Vector4[] in_homoClipPlaneNormals, float[,] in_dotClipPlanesArray, int [] in_clipCodeArray) {
        ClipStepOutput output;

        List<Vector4> out_hxhyhzwList = new List<Vector4>();
        List<Vector2> out_uvList = new List<Vector2>();
        List<int> out_triListIndexList = new List<int>();

        // Process one triangle at a time
        int triCount = in_triListIndexArray.Length / 3;
        for (int iTri = 0 ; iTri < triCount; ++iTri) {
            // Determine if we can trivially accept or reject the triangle based on the clip code
            int clipCodeAND = 0;
            int clipCodeOR = 0;
            for (int i = 0; i < 3; ++i) {
                int index = in_triListIndexArray[3 * iTri + i];
                int clipCode = in_clipCodeArray[index];
                if (i == 0) clipCodeAND = clipCode; else clipCodeAND &= clipCode;
                clipCodeOR |= clipCode;
            }
            if (clipCodeAND != 0) {
                // All vertices are outside at least one of the planes, we can trivially reject this triangle
            } else if (clipCodeOR == 0) {
                // All vertices are inside all the planes, we can trivially accept this triangle
                int baseIndex = out_hxhyhzwList.Count;
                for (int i = 0; i < 3; ++i) {
                    int index = in_triListIndexArray[3 * iTri + i];
                    out_hxhyhzwList.Add(in_hxhyhzwArray[index]);
                    out_uvList.Add(in_uvArray[index]);
                    out_triListIndexList.Add(baseIndex + i);
                }
            } else {
                // Triangle needs to be clipped.  The result of the clip will be a convex polygon.
                List<Vector4> polygonHomoPosList = new List<Vector4>();
                List<Vector2> polygonUVList = new List<Vector2>();
                for (int i = 0; i < 3; ++i) {
                    int index = in_triListIndexArray[3 * iTri + i];
                    polygonHomoPosList.Add(in_hxhyhzwArray[index]);
                    polygonUVList.Add(in_uvArray[index]);
                }

                // Consider one clip plane at a time
                for (int iPlane = 0; iPlane < in_homoClipPlaneNormals.Length; ++iPlane) {
                    // We can skip any planes that did not clip the original triangle
                    if ((clipCodeOR & (1 << iPlane)) != 0) {
                        // Clip the polygon and replace our polygon data with the new data
                        ClipPolygonOutput clipPolygonOutput = ClipPolygon(polygonHomoPosList, polygonUVList, in_homoClipPlaneNormals[iPlane]);
                        polygonHomoPosList = clipPolygonOutput.hxhyhzwList;
                        polygonUVList = clipPolygonOutput.uvList;
                    }
                }

                // Add clipped geometry to output and create a triangle fan for the polygon
                int baseIndex = out_hxhyhzwList.Count;
                out_hxhyhzwList.AddRange(polygonHomoPosList);
                out_uvList.AddRange(polygonUVList);
                for (int i = 0; i < polygonHomoPosList.Count - 2; ++i) {
                    out_triListIndexList.Add(baseIndex);
                    out_triListIndexList.Add(baseIndex + i + 1);
                    out_triListIndexList.Add(baseIndex + i + 2);
                }

            }
        }
    
        // Convert our lists to array
        output.hxhyhzwArray = out_hxhyhzwList.ToArray();
        output.uvArray = out_uvList.ToArray();
        output.triListIndexArray = out_triListIndexList.ToArray();
        return output;
    }

    struct RasterSetupStepOutput {
        // Screen x, screen y and reciprocal of w
        public Vector3[] sxsyrhwArray;
        // Reciprocal of w times uv
        public Vector2[] rhwuvArray;
    }

    // Compute rasterization interpolants
    static RasterSetupStepOutput RasterSetupStep(Vector2 screenSizeInPixels, Vector4[] hxhyhzwArray, Vector2[] uvArray) {

        // Alloc output
        RasterSetupStepOutput output;
        int count = hxhyhzwArray.Length;
        output.sxsyrhwArray = new Vector3[count];
        output.rhwuvArray = new Vector2[count];

        for (int i = 0; i < count; ++i) {
            // Compute screen x,y and rhw
            Vector4 hxhyhzw = hxhyhzwArray[i];
            Vector4 sxsyrhw = output.sxsyrhwArray[i];
            sxsyrhw.x = (hxhyhzw.x/hxhyhzw.w + 1f)*.5f * screenSizeInPixels.x;
            sxsyrhw.y = (hxhyhzw.y/hxhyhzw.w + 1f)*.5f * screenSizeInPixels.y;
            sxsyrhw.z = 1f / hxhyhzw.w;

            // Compute reciprocal w times uv
            Vector2 uv = uvArray[i];
            Vector2 rhwuv = output.rhwuvArray[i];
            rhwuv = sxsyrhw.z * uv;
        }

        return output;
    }

    // Barycentric coords for the input position given three points of the triangle
    static float[] ComputeBarycentricCoords(float px, float py, float v1x, float v1y, float v2x, float v2y, float v3x, float v3y)
    {
        float w1 = ((v2y - v3y) * (px - v3x) + (v3x - v2x) * (py - v3y)) /  ((v2y - v3y) * (v1x - v3x) + (v3x - v2x) * (v1y - v3y));
        float w2 = ((v3y - v1y) * (px - v3x) + (v1x - v3x) * (py - v3y)) /  ((v2y - v3y) * (v1x - v3x) + (v3x - v2x) * (v1y - v3y));
        float w3 = 1.0f - w1 - w2;
        return new float[] {w1, w2, w3};
    }

    // Interpolate the UV coordinate for the given point using barycentric coords
    static Vector2 ComputeInterpolatedUV(float sx, float sy, Vector3[] sxsyrhwArray, Vector2[] rhwuvArray) {
        float[] w = ComputeBarycentricCoords(sx, sy, sxsyrhwArray[0].x, sxsyrhwArray[0].y, sxsyrhwArray[1].x, sxsyrhwArray[1].y, sxsyrhwArray[2].x, sxsyrhwArray[2].y);
        float rhw = w[0] * sxsyrhwArray[0].z + w[1] * sxsyrhwArray[1].z + w[2] * sxsyrhwArray[2].z;
        Vector2 rhwuv = w[0] * rhwuvArray[0] + w[1] * rhwuvArray[1] + w[2] * rhwuvArray[2];
        Vector2 uv = rhwuv / rhw;
        return uv;
    }

    struct ComputeMipLevelsOutput {
        public float[] mipLevelArray;
    }

    // Compute mip levels for every vertex of every triangle.
    // There really should only be one mip level per vertex, but because of the way the geometry is specified, we can't guarantee this
    static ComputeMipLevelsOutput ComputeMipLevels(Vector2 in_textureSizeInPixels, Vector3[] in_sxsyrhwArray, Vector2[] in_rhwuvArray, int[] in_triListIndexArray) {
        // Allocate output
        ComputeMipLevelsOutput output;
        output.mipLevelArray = new float[in_triListIndexArray.Length];

        // Process one triangle at a time
        int count = in_triListIndexArray.Length / 3;
        for (int iTri = 0; iTri < count; ++iTri) {
            // Retrieve triangle info for easier processing
            Vector3[] sxsyrhwArray = new Vector3[3];
            Vector2[] rhwuvArray = new Vector2[3];
            for (int i = 0; i < 3; ++i) { 
                int index = in_triListIndexArray[iTri * 3 + i];
                sxsyrhwArray[i] = in_sxsyrhwArray[index];
                rhwuvArray[i] = in_rhwuvArray[index];
            }

            // Process each point of the triangle
            const float d = 1.0f; // Distance from each point to interpolate
            for (int i = 0; i < 3; ++i) {
                // Data for this point
                float sx = sxsyrhwArray[i].x;
                float sy = sxsyrhwArray[i].y;

                // Compute interpolated UV for d pixel in both screen directions
                Vector2 uv_dx = ComputeInterpolatedUV(sx + d, sy, sxsyrhwArray, rhwuvArray);
                Vector2 uv_dy = ComputeInterpolatedUV(sx, sy + d, sxsyrhwArray, rhwuvArray);

                // Complete computation of partial derivatives
                Vector2 uv = rhwuvArray[i] / sxsyrhwArray[i].z;
                Vector2 duv_dx = uv_dx - uv;
                Vector2 duv_dy = uv_dy - uv;

                // Compute uv area of parallelogram mapped by a single pixel, multiply by texture dim, sqrt for mip level! finally!                
                float a_uv = Mathf.Abs(duv_dx.x * duv_dy.y - duv_dy.x * duv_dx.y) / (d * d);
                float mipLevel = Mathf.Sqrt(a_uv * in_textureSizeInPixels.x * in_textureSizeInPixels.y);
                output.mipLevelArray[iTri * 3 + i] = mipLevel;
            }
        }
        return output;
    }

    // Normals of clip planes in homogeneous space facing inwards for LH space
    static Vector4[] homoClipPlaneNormals = new Vector4[] {
        new Vector4( 1.0f, 0.0f, 0.0f, 1.0f),  // x left
        new Vector4(-1.0f, 0.0f, 0.0f, 1.0f),  // x right
        new Vector4( 0.0f,-1.0f, 0.0f, 1.0f),  // y top
        new Vector4( 0.0f, 1.0f, 0.0f, 1.0f),  // y bottom
        new Vector4( 0.0f, 0.0f, 1.0f, 1.0f),  // z near
        new Vector4( 0.0f, 0.0f,-1.0f, 1.0f),  // z far
    };

    public struct DoOutput {
        public float minMipLevel;
        public float maxMipLevel;
        public bool isOutsideView;
    };

    public static DoOutput Do(Matrix4x4 in_worldViewProjMatrix, Vector2 in_screenSizeInPixels, Vector2 in_textureSizeInPixels, Vector3[] in_xyzArray, Vector2[] in_uvArray, int[] in_triListIndexArray) {
        DoOutput output;
        output.isOutsideView = false;
        output.minMipLevel = float.MaxValue;
        output.maxMipLevel = float.MinValue;

        // Transform all mesh into homogeneous space and get clip codes
        TransformStepOutput transformStepOutput = TransformStep(in_worldViewProjMatrix, homoClipPlaneNormals, in_xyzArray);

        // Trivially accept or reject based on clip AND and OR
        if (transformStepOutput.clipCodeAND != 0) {
            // The entire set of vertices is outside the view, so we can early exit
            output.isOutsideView = true;
        } else {
            ClipStepOutput clipStepOutput;
            if (transformStepOutput.clipCodeOR == 0) {
                // trivially accept
                clipStepOutput.hxhyhzwArray = transformStepOutput.hxhyhzwArray;
                clipStepOutput.uvArray = in_uvArray;
                clipStepOutput.triListIndexArray = in_triListIndexArray;
            } else {
                // If we cannot trivially accept or reject, clip
                clipStepOutput = ClipStep(transformStepOutput.hxhyhzwArray, in_uvArray, in_triListIndexArray, homoClipPlaneNormals, transformStepOutput.dotClipPlanesArray, transformStepOutput.clipCodeArray);
            }

            // Compute interpolants we will use for our mip calculation
            RasterSetupStepOutput rasterSetupStepOutput = RasterSetupStep(in_screenSizeInPixels, clipStepOutput.hxhyhzwArray, clipStepOutput.uvArray);

            // Compute mip levels for every vertex of every triangle
            ComputeMipLevelsOutput computeMipLevelsOutput = ComputeMipLevels(in_textureSizeInPixels, rasterSetupStepOutput.sxsyrhwArray, rasterSetupStepOutput.rhwuvArray, clipStepOutput.triListIndexArray);

            // Compute max/min values for return
            for (int i = 0; i < computeMipLevelsOutput.mipLevelArray.Length; ++i) {
                float x = computeMipLevelsOutput.mipLevelArray[i];
                output.maxMipLevel = Mathf.Max(output.maxMipLevel, x);
                output.minMipLevel = Mathf.Min(output.minMipLevel, x);
            }
        }
        return output;
    }
}
