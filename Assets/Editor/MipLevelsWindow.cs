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

    // Sample for a single texture
    struct TextureSample {
        public Texture texture;
        // min and max mip levels used on visible regions of geometry using this texture 
        public float minMipLevel;
        public float maxMipLevel;
        // true if all uses of this texture are outside the view
        public bool isOutsideView;
        // Sorted array of sub mesh samples.  Sorted from highest min value to lowest
        public SubMeshSample[] subMeshSampleArray;
    }

    // Sample data for each sub mesh
    struct SubMeshSample {
        // Object, mesh and sub mesh source
        public GameObject gameObject;
        public int subMeshIndex;
        public Mesh mesh;
        public MeshRenderer meshRenderer;
        // Mip level data
        public ComputeMipLevels.MainOutput computeMipLevelsMainOutput;
    }

    // Sorted array of samples for each texture.  Sorted from highest min values to lowest
    TextureSample[] textureSampleArray;

    // Analyze the scene and save output to textureSampleArray
    void Analyze()
    {
        // Use the main camera
        Camera camera = Camera.main;
        Vector2 screenSizeInPixels = new Vector2(camera.pixelWidth, camera.pixelHeight);
        Matrix4x4 viewProjMatrix = camera.projectionMatrix * camera.worldToCameraMatrix;

        // Build our output here
        Dictionary<Texture, List<SubMeshSample>> textureToSubMeshSampleList = new Dictionary<Texture, List<SubMeshSample>>();

        // Walk all game objects that are active in the scene
        GameObject[] allObjects = Object.FindObjectsOfType<GameObject>() ;
        foreach(GameObject gameObject in allObjects) {
            if (gameObject.activeInHierarchy) {
                // Look for a mesh 
                MeshFilter meshFilter = gameObject.GetComponent<MeshFilter>();
                Mesh mesh = (meshFilter == null) ? null : meshFilter.sharedMesh;
                MeshRenderer meshRenderer = gameObject.GetComponent<MeshRenderer>();
                if (mesh != null && meshRenderer != null) {
                    // Walk all submeshes
                    Matrix4x4 worldViewProjMatrix = viewProjMatrix * meshRenderer.localToWorldMatrix;
                    for (int subMeshIndex = 0; subMeshIndex < mesh.subMeshCount; ++subMeshIndex) {
                        // Look for materials that use textures
                        Material material = meshRenderer.sharedMaterials[subMeshIndex];
                        Texture texture = material.mainTexture;
                        if (texture != null) {
                            // Find or create data for this texture
                            List<SubMeshSample> subMeshSampleList = null;
                            if (textureToSubMeshSampleList.TryGetValue(texture, out subMeshSampleList) == false) {
                                subMeshSampleList = new List<SubMeshSample>();
                                textureToSubMeshSampleList.Add(texture, subMeshSampleList);
                            }
                            // Compute the mip levels for this submesh
                            int[] triangles = mesh.GetTriangles(subMeshIndex);
                            Vector2 textureSizeInPixels = new Vector2(texture.width, texture.height);
                            ComputeMipLevels.MainOutput computeMipLevelsMainOutput = ComputeMipLevels.Main(worldViewProjMatrix, screenSizeInPixels, textureSizeInPixels, mesh.vertices, mesh.uv, triangles);                        
                            SubMeshSample subMeshSample;
                            subMeshSample.gameObject = gameObject;
                            subMeshSample.subMeshIndex = subMeshIndex;
                            subMeshSample.mesh = mesh;
                            subMeshSample.meshRenderer = meshRenderer;
                            subMeshSample.computeMipLevelsMainOutput = computeMipLevelsMainOutput;
                            // Record
                            subMeshSampleList.Add(subMeshSample);
                        }
                    }
                }
            }
        }

        // Convert our build data to our final representation
        int numTextureSample = textureToSubMeshSampleList.Count;
        textureSampleArray = new TextureSample[numTextureSample];
        int iTextureSample = 0;
        foreach (KeyValuePair<Texture, List<SubMeshSample> > keyValuePair in textureToSubMeshSampleList) {
            Texture texture = keyValuePair.Key;
            List<SubMeshSample> subMeshSampleList = keyValuePair.Value;
            int numSubMeshSample = subMeshSampleList.Count;
            TextureSample textureSample;
            textureSample.texture = texture;
            textureSample.minMipLevel = float.MaxValue;
            textureSample.maxMipLevel = float.MinValue;
            textureSample.isOutsideView = true;            
            textureSample.subMeshSampleArray = new SubMeshSample[numSubMeshSample];
            int iSubMesh = 0;
            foreach(SubMeshSample subMeshSample in subMeshSampleList) {
                textureSample.minMipLevel = Mathf.Min(textureSample.minMipLevel, subMeshSample.computeMipLevelsMainOutput.minMipLevel);
                textureSample.maxMipLevel = Mathf.Max(textureSample.maxMipLevel, subMeshSample.computeMipLevelsMainOutput.maxMipLevel);
                textureSample.isOutsideView &= subMeshSample.computeMipLevelsMainOutput.isOutsideView;
                textureSample.subMeshSampleArray[iSubMesh] = subMeshSample;
                ++iSubMesh;
            }
            textureSampleArray[iTextureSample] = textureSample;
            ++iTextureSample;
        }

        // Sort all the sub mesh samples by min mip level (also isOutsideView)
        for (int i = 0; i < textureSampleArray.Length; ++i) {
            System.Array.Sort<SubMeshSample>(textureSampleArray[i].subMeshSampleArray, delegate(SubMeshSample lhs, SubMeshSample rhs) {
                if (!lhs.computeMipLevelsMainOutput.isOutsideView &&  rhs.computeMipLevelsMainOutput.isOutsideView) return -1;
                if ( lhs.computeMipLevelsMainOutput.isOutsideView && !rhs.computeMipLevelsMainOutput.isOutsideView) return  1;
                if (lhs.computeMipLevelsMainOutput.minMipLevel < rhs.computeMipLevelsMainOutput.minMipLevel) return -1;
                if (lhs.computeMipLevelsMainOutput.minMipLevel > rhs.computeMipLevelsMainOutput.minMipLevel) return 1;
                return 0;
            });
        }
        // Sort all the texture samples by min mip level with largest first (also isOutSideView)
        System.Array.Sort<TextureSample>(textureSampleArray, delegate(TextureSample lhs, TextureSample rhs) {
            if (!lhs.isOutsideView &&  rhs.isOutsideView) return -1;
            if ( lhs.isOutsideView && !rhs.isOutsideView) return  1;
            if (lhs.minMipLevel < rhs.minMipLevel) return 1;
            if (lhs.minMipLevel > rhs.minMipLevel) return -1;
            return 0;
        });
    }

    // The previous GUI item values
    Vector2 prevScrollPosTextures;
    int prevSelGridIndexTexture = -1;
    Vector2 prevScrollPosSubMesh;
    int prevSelGridIndexSubmesh = -1;
 
    // Draw our GUI and process user input
    void OnGUI() {
        // Button to analyze
        if (GUILayout.Button("Analyze", GUILayout.ExpandWidth(false))) {
            Analyze();
            // Reset our UI
            prevScrollPosTextures = new Vector2(0.0f, 0.0f);
            prevSelGridIndexTexture = -1;
            prevScrollPosSubMesh = new Vector2(0.0f, 0.0f);
            prevSelGridIndexSubmesh = -1;
                    
        }
        GUILayout.Space(10);

        if (textureSampleArray != null) {

            // Scroll view of textures
            int heightScrollViewTexture = (int)(0.3f * position.height);
            prevScrollPosTextures = EditorGUILayout.BeginScrollView(prevScrollPosTextures, GUILayout.Height(heightScrollViewTexture));
            GUIContent[] guiContentArrayTexture = new GUIContent[textureSampleArray.Length];
            for (int i = 0; i < textureSampleArray.Length; ++i) {
                TextureSample textureSample = textureSampleArray[i];                
                GUIContent guiContent = new GUIContent("T \"" + textureSample.texture.name + "\" -> " +
                    (!textureSample.isOutsideView ? "min mip: " + System.String.Format("{0:0.#}", textureSample.minMipLevel) + " #SM " + textureSample.subMeshSampleArray.Length : "OUTSIDE"));
                guiContentArrayTexture[i] = guiContent;
            }
            int curSelGridIndexTexture = GUILayout.SelectionGrid(prevSelGridIndexTexture, guiContentArrayTexture, 1);
            if (curSelGridIndexTexture != prevSelGridIndexTexture) {
                // Selected texture has changed
                Selection.activeObject = textureSampleArray[curSelGridIndexTexture].texture;                
                prevSelGridIndexTexture = curSelGridIndexTexture;
                prevSelGridIndexSubmesh = -1;
            }
            EditorGUILayout.EndScrollView();

            GUILayout.Space(10);

            if (prevSelGridIndexTexture >= 0) {
                // Scroll view of sub meshes for currently selected texture
                int heightScrollViewSubmesh = (int)(0.3f * position.height);
                prevScrollPosSubMesh = EditorGUILayout.BeginScrollView(prevScrollPosSubMesh, GUILayout.Height(heightScrollViewSubmesh));
                SubMeshSample[] subMeshSampleArray = textureSampleArray[prevSelGridIndexTexture].subMeshSampleArray;
                GUIContent[] guiContentArraySubMesh = new GUIContent[subMeshSampleArray.Length];
                for (int i = 0; i < subMeshSampleArray.Length; ++i) {
                    SubMeshSample subMeshSample = subMeshSampleArray[i];
                    GUIContent guiContent = new GUIContent("GO \"" + subMeshSample.gameObject.name + "\" / M \"" + subMeshSample.mesh.name + "\" / SM " + subMeshSample.subMeshIndex +
                        (!subMeshSample.computeMipLevelsMainOutput.isOutsideView ? 
                            " -> min mip: " + System.String.Format("{0:0.#}", subMeshSample.computeMipLevelsMainOutput.minMipLevel) : "OUTSIDE"));
                    guiContentArraySubMesh[i] = guiContent;                
                }
                int curSelGridIndexSubMesh = GUILayout.SelectionGrid(prevSelGridIndexSubmesh, guiContentArraySubMesh, 1);
                if (curSelGridIndexSubMesh != prevSelGridIndexSubmesh) {
                    // Selected submesh has changed
                    Selection.activeTransform = subMeshSampleArray[curSelGridIndexSubMesh].gameObject.transform;
                    prevSelGridIndexSubmesh = curSelGridIndexSubMesh;
                }
                EditorGUILayout.EndScrollView();
            }
        }
    }

}
