using System;
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using Unity.VisualScripting;
using Unity.VisualScripting.FullSerializer;
using UnityEngine;

//TGParser stands for Tetgen Parser. Yes, cringe name, problem, officer?
//GNU stands for GNU is Not Linux and nobody makes a big deal out of it.
//If this were a bad naming competition, I'm sure I wouldn't even make it to top 10.

public static class TGParser
{
    private static readonly CultureInfo locale = new CultureInfo("en-US");
    private static readonly string[] split_separators_values = {"\n", "\t", "\r", "\t", " ", "\v", "\r\n"};
    private static readonly string[] split_separators_lines = { "\n", "\r", "\r\n" };
    private static readonly StringSplitOptions split_options = StringSplitOptions.RemoveEmptyEntries;

    public static bool flipX = true;
    public static bool flipY = false;
    public static bool flipZ = false;

    private static string[] FilterComments(string text)
    {
        // generate the lines and remove lines that start with the symbol "#"
        string[] lines = text.Split(split_separators_lines, split_options).Where(line => !line.TrimStart().StartsWith("#")).ToArray();

        // remove the contents of each line from the moment they start to contain a "#" symbol to support inline comments
        for(int k = 0; k < lines.Length; ++k)
        {
            for (int i = 0; i < lines[k].Length; ++i)
            {
                if (lines[k][i] == '#')
                {
                    lines[k] = lines[k].Substring(0, i); // generates a sub string starting at index 0 and with length i (since index i has length i + 1, then the generated string is of the expected length when removing everything after the comment)
                }
            }
        }
        
        return lines;
    }

    private static string[] FileToLines(TextAsset file)
    {
        string[] lines = FilterComments(file.text);
        return lines;
    }

    private static string[] LineToValues(string line)
    {
        string[] values = line.Split(split_separators_values, split_options);
        return values;
    }

    public static float GetValueFloat(string s)
    {
        float val = float.Parse(s, TGParser.locale);
        return val;
    }
    public static int GetValueInt(string s)
    {
        int val = int.Parse(s, TGParser.locale);
        return val;
    }

    public static Node[] ParseFileNode(TextAsset file_node, Transform transform)
    {
        string[] lines = FileToLines(file_node);
        int num_nodes = ReadLength(lines[0]);
        Node[] nodes = new Node[num_nodes];

        for (int i = 0; i < num_nodes; ++i)
        {
            Vector3 current_position = transform.TransformPoint(ReadPositionFlipped(lines[i + 1]));
            Node current_node = new Node(current_position);
            nodes[i] = current_node;
            nodes[i].mass = 0.0f; // Set the starting mass to 0 so that we can properly add it through the mass addition loop in the elastic solid class.
        }

        return nodes;
    }

    public static ElasticSolid.Edge[] ParseFileEdge(TextAsset file_edges, Transform transform)
    {
        string[] lines = FileToLines(file_edges);
        int num_edges = ReadLength(lines[0]);
        ElasticSolid.Edge[] edges = new ElasticSolid.Edge[num_edges];

        for (int i = 0; i < num_edges; ++i)
        {
            ElasticSolid.Edge edge = ReadEdge(lines[i + 1]);
            edges[i] = edge;
        }

        return edges;
    }

    public static ElasticSolid.Tetrahedron[] ParseFileEle(TextAsset file_ele, Transform transform, Node[] nodes)
    {
        string[] lines = FileToLines(file_ele);
        int num_tetrahedra = ReadLength(lines[0]);
        ElasticSolid.Tetrahedron[] tetrahedra = new ElasticSolid.Tetrahedron[num_tetrahedra];

        for (int i = 0; i < num_tetrahedra; ++i)
        {
            ElasticSolid.Tetrahedron current_tetrahedron = ReadTetrahedron(lines[i + 1], nodes);
            tetrahedra[i] = current_tetrahedron;
        }

        return tetrahedra;
    }

    private static int ReadLength(string line)
    {
        string[] values = LineToValues(line);
        int ans = GetValueInt(values[0]);
        return ans;
    }

    private static Vector3 ReadPosition(string line)
    {
        string[] values = LineToValues(line);
        float x, y, z;
        x = GetValueFloat(values[1]);
        y = GetValueFloat(values[2]);
        z = GetValueFloat(values[3]);
        Vector3 ans = new Vector3(x,y,z);
        return ans;
    }

    private static Vector3 ReadPositionFlipped(String line)
    {
        Vector3 pos = ReadPosition(line);
        if (flipX) pos.x *= -1;
        if (flipY) pos.y *= -1;
        if (flipZ) pos.z *= -1;
        return pos;
    }

    private static ElasticSolid.Edge ReadEdge(string line)
    {
        string[] values = LineToValues(line);
        int a = GetValueInt(values[1]) - 1;
        int b = GetValueInt(values[2]) - 1;
        ElasticSolid.Edge ans = new ElasticSolid.Edge(a, b);
        return ans;
    }

    private static ElasticSolid.Tetrahedron ReadTetrahedron(string line, Node[] nodes)
    {
        string[] values = LineToValues(line);
        int node1 = GetValueInt(values[1]) - 1;
        int node2 = GetValueInt(values[2]) - 1;
        int node3 = GetValueInt(values[3]) - 1;
        int node4 = GetValueInt(values[4]) - 1;
        ElasticSolid.Tetrahedron ans = new ElasticSolid.Tetrahedron(node1, node2, node3, node4, nodes);
        return ans;
    }

}
