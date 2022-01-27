using UnityEngine;
using UnityEditor;
using System.Collections;
using System.Collections.Generic;

public class MipLevelsWindow : EditorWindow 
{
    [MenuItem( "Window/Mip Levels Window" )]
    static void OpenWindow()
    {
        EditorWindow.GetWindow<MipLevelsWindow>( false, "Mip Levels" );
    }

    void OnEnable()
    {
        EditorApplication.update -= OnUpdate;
        EditorApplication.update += OnUpdate;
    }

    void OnDisable()
    {
        EditorApplication.update -= OnUpdate;
    }

    void OnUpdate()
    {
    }

    struct SubMeshSample {
        public int subMeshIndex;
        public Mesh mesh;
        public MeshRenderer meshRenderer;
        public ComputeMipLevels.MainOutput computeMipLevelsMainOutput;
    }

    SubMeshSample[] subMeshSampleArray;

    void Sample()
    {
        Camera camera = Camera.main;
        Vector2 screenSizeInPixels = new Vector2(camera.pixelWidth, camera.pixelHeight);
        Matrix4x4 viewProjMatrix = camera.projectionMatrix * camera.worldToCameraMatrix;

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
        MeshFilter[] meshFilterArray = listMeshFilter.ToArray();
        MeshRenderer[] meshRendererArray = listMeshRenderer.ToArray();

        List<SubMeshSample> subMeshSampleList = new List<SubMeshSample>();
        for (int meshIndex = 0; meshIndex < meshFilterArray.Length; ++meshIndex) {
            MeshFilter meshFilter = meshFilterArray[meshIndex];
            MeshRenderer meshRenderer = meshRendererArray[meshIndex];
            Mesh mesh = meshFilter.mesh;
            Matrix4x4 worldViewProjMatrix = viewProjMatrix * meshRenderer.localToWorldMatrix;
            Material[] materials = meshRenderer.materials;
            if (mesh.subMeshCount <= materials.Length) {
                for (int subMeshIndex = 0; subMeshIndex < mesh.subMeshCount; ++subMeshIndex) {
                    Material material = materials[subMeshIndex];
                    int[] triangles = mesh.GetTriangles(subMeshIndex);
                    if (material.mainTexture != null) {
                        Vector2 textureSizeInPixels = new Vector2(material.mainTexture.width, material.mainTexture.height);
                        ComputeMipLevels.MainOutput computeMipLevelsMainOutput = ComputeMipLevels.Main(worldViewProjMatrix, screenSizeInPixels, textureSizeInPixels, mesh.vertices, mesh.uv, triangles);                        
                        SubMeshSample subMeshSample;
                        subMeshSample.subMeshIndex = subMeshIndex;
                        subMeshSample.mesh = mesh;
                        subMeshSample.meshRenderer = meshRenderer;
                        subMeshSample.computeMipLevelsMainOutput = computeMipLevelsMainOutput;
                        subMeshSampleList.Add(subMeshSample);
                    }
                }
            }
        }
        subMeshSampleArray = subMeshSampleList.ToArray();
    }

    Vector2 prevScrollPos;
    int prevSelGridIndex = 0;
 
    void OnGUI() {

        if (GUILayout.Button("Sample", GUILayout.ExpandWidth(false))) {
            Sample();
        }
        GUILayout.Space(10);
        prevScrollPos = EditorGUILayout.BeginScrollView(prevScrollPos); //, GUILayout.Width(100), GUILayout.Height(100));
        if (subMeshSampleArray != null) {
            GUIContent[] guiContentArray = new GUIContent[subMeshSampleArray.Length];
            for (int i = 0; i < subMeshSampleArray.Length; ++i) {
                SubMeshSample subMeshSample = subMeshSampleArray[i];
                Texture texture = subMeshSample.meshRenderer.materials[subMeshSample.subMeshIndex].mainTexture;
                GUIContent guiContent = new GUIContent("M / " + subMeshSample.mesh.name + "/ SM " + subMeshSample.subMeshIndex + "/ T " + texture.name + " -> min/max: " + 
                    System.String.Format("{0:0.#}", subMeshSample.computeMipLevelsMainOutput.minMipLevel) + "/" + 
                    System.String.Format("{0:0.#}", subMeshSample.computeMipLevelsMainOutput.maxMipLevel)  );
                guiContentArray[i] = guiContent;
            }
            int curSelGridIndex = GUILayout.SelectionGrid(prevSelGridIndex, guiContentArray, 1);
            if (curSelGridIndex != prevSelGridIndex) {
                Debug.Log(prevSelGridIndex + " to " + curSelGridIndex);
                prevSelGridIndex = curSelGridIndex;
            }
        }
        EditorGUILayout.EndScrollView();

    }

}
