using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Drawing;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Linq;
using MathNet.Numerics.LinearAlgebra.Factorization;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using MathNet.Numerics.Integration;



namespace Master.Components
{
    public class Assembly3DBeamLinear : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Assembly2DTruss class.
        /// </summary>
        public Assembly3DBeamLinear()
          : base("Assembly3DBeamLinear", "Nickname",
              "Description",
              "Panda", "3DBeam")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Beams", "B", "BeamClass object", GH_ParamAccess.list);
            pManager.AddGenericParameter("Boundary Conditions", "BC", "BcClass object", GH_ParamAccess.list);
            pManager.AddGenericParameter("Loads", "L", "LoadClass object", GH_ParamAccess.list);
            pManager.AddNumberParameter("x-position", "Pos", "x-Position", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Displacement [mm] in x", "def", "Displacement [mm]", GH_ParamAccess.list);
            pManager.AddPointParameter("Rotation [rad] in x", "R", "Forces in points", GH_ParamAccess.list);
            pManager.AddPointParameter("Forces [N]", "F", "Forces in points", GH_ParamAccess.list);
            pManager.AddPointParameter("Moment [Nmm]", "R", "Rotation in points", GH_ParamAccess.list);
            pManager.AddPointParameter("DisplacementOfNodes", "displ", "Deformed geometry", GH_ParamAccess.list);
            pManager.AddLineParameter("Deformed Geometry", "defGeo", "Deformed geometry", GH_ParamAccess.list);
            pManager.AddNumberParameter("Strain in x", "E", "strain ", GH_ParamAccess.list);
            pManager.AddNumberParameter("Stress [N/mm^2] in x", "S", "stress [N/mm^2] ", GH_ParamAccess.list);
            pManager.AddNumberParameter("Maximum moment", "M_max", "Maximum moment ", GH_ParamAccess.item);
            pManager.AddNumberParameter("umax", "Umax", "maximum displacement ", GH_ParamAccess.item);
            pManager.AddLineParameter("moment diagram", "ben", "to see the moment ", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {


            //input


            List<BeamClassLinear> bars = new List<BeamClassLinear>();
            List<LoadClass> lc = new List<LoadClass>();
            List<BcClass> bcc = new List<BcClass>();
            double x = new double();


            DA.GetDataList(0, bars);
            DA.GetDataList(1, bcc);
            DA.GetDataList(2, lc);
            DA.GetData(3, ref x);


            List<Point3d> pts = new List<Point3d>();
            List<NodeClass> nodes = new List<NodeClass>();
            foreach (BeamClassLinear b in bars)
            {

                if (!pts.Contains(b.startNode.pt))
                {
                    pts.Add(b.startNode.pt);
                    nodes.Add(b.startNode);
                }
                if (!pts.Contains(b.endNode.pt))
                {
                    pts.Add(b.endNode.pt);
                    nodes.Add(b.endNode);
                }
            }

            List<int> BCList = CreateBCList(bcc, nodes);

            Vector<double> LoadList = CreateLoadList(lc, nodes);



            CreateGlobalStiffnesMatrix(bars, pts, out Matrix<double> k_tot, out List<Matrix<double>> Lk_eg, out List<Matrix<double>> T_list, out List<Matrix<double>> T_list_six);

            var R = Vector<double>.Build.DenseOfArray(LoadList.ToArray());

            CreateReducedStiffnesMatrix(BCList, k_tot, R, out Matrix<double> K_red);

            Matrix<double> invK = K_red.Inverse();

            var def = invK.Multiply(R);
            var displNodes = new List<Point3d>();

            //Vector<double> def = K_red.Cholesky().Solve(R);     //r


            CreateGenerelizedShapeFunc(bars, k_tot, x, out Matrix<double> N, out Matrix<double> dN);


            CreateForces(bars, pts, def, Lk_eg, N, dN, T_list, T_list_six, out Vector<double> forces, out Vector<double> rotation, out Vector<double> strain, out Vector<double> stress, out List<Vector<double>> u, out double M_max, out double umax, out List<Line> crv);


            for (int i = 0; i < nodes.Count; i++)
            {
                Point3d node = nodes[i].pt;
                int id = nodes[i].Id;
                displNodes.Add(new Point3d(node.X + def[6 * id] / 1000, node.Y + def[6 * id + 1] / 1000, node.Z + def[6 * id + 2] / 1000));
            }


            //lage lister med displacement

            List<Point3d> force_lst = new List<Point3d>();
            List<Point3d> mom_lst = new List<Point3d>();
            List<Point3d> disp_lst = new List<Point3d>();
            List<Point3d> rot_lst = new List<Point3d>();

            //gjør om små tall til 0
            for (int i = 0; i < forces.Count; i++)
            {
                if (Math.Abs(forces[i]) <= 0.00001)
                {
                    forces[i] = 0;
                }


                if (Math.Abs(rotation[i]) <= 0.00001)
                {
                    rotation[i] = 0;
                }

            }

            foreach (Vector uu in u)
            {
                disp_lst.Add(new Point3d(uu[0], uu[1], uu[2]));
                rot_lst.Add(new Point3d(uu[3], uu[4], uu[5]));
            }

            // endrer output til liste med punkter
            for (int i = 0; i < nodes.Count; i++)
            {


                if (BCList.Contains(nodes[i].Id * 6) | BCList.Contains(nodes[i].Id * 6 + 1) | BCList.Contains(nodes[i].Id * 6 + 2))
                {
                    force_lst.Add(new Point3d(forces[nodes[i].Id * 3] - LoadList[nodes[i].Id * 6], forces[nodes[i].Id * 3 + 1] - LoadList[nodes[i].Id * 6 + 1], forces[nodes[i].Id * 3 + 2] - LoadList[nodes[i].Id * 6 + 2]));
                    mom_lst.Add(new Point3d(rotation[nodes[i].Id * 3] + LoadList[nodes[i].Id * 6 + 3], rotation[nodes[i].Id * 3 + 1] + LoadList[nodes[i].Id * 6 + 4], rotation[nodes[i].Id * 3 + 2] + LoadList[nodes[i].Id * 6 + 5]));

                }


            }

            //Lager deformed geometry

            List<Line> deformedgeometry = new List<Line>();

            foreach (BeamClassLinear b in bars)
            {
                Point3d p1 = new Point3d(b.startNode.pt.X + (def[b.startNode.Id * 6] / 1000.00), b.startNode.pt.Y + (def[b.startNode.Id * 6 + 1] / 1000.00), b.startNode.pt.Z + (def[b.startNode.Id * 6 + 2] / 1000.00));
                Point3d p2 = new Point3d(b.endNode.pt.X + (def[b.endNode.Id * 6] / 1000.00), b.endNode.pt.Y + (def[b.endNode.Id * 6 + 1] / 1000.00), b.endNode.pt.Z + (def[b.endNode.Id * 6 + 2] / 1000.00));
                Line line1 = new Line(p1, p2);
                

                deformedgeometry.Add(line1);
                

            }



            //output
            DA.SetDataList(0, disp_lst);
            DA.SetDataList(1, rot_lst);
            DA.SetDataList(2, force_lst);
            DA.SetDataList(3, mom_lst);
            DA.SetDataList(4, displNodes);
            DA.SetDataList(5, deformedgeometry);
            DA.SetDataList(6, strain);
            DA.SetDataList(7, stress);
            DA.SetData(8, M_max);
            DA.SetData(9, umax);
            DA.SetDataList(10, crv);



        }

        private List<int> CreateBCList(List<BcClass> _BcValue, List<NodeClass> _nodes)  //making a list with indexes of fixed BC
        {

            //List<int> BCs = new List<int>();
            List<int> BCsIndex = new List<int>();


            for (int i = 0; i < _nodes.Count; i++)
            {
                Point3d node = _nodes[i].pt;
                int id = _nodes[i].Id;
                foreach (var b in _BcValue)
                {
                    if (b.Coordinate.DistanceTo(node) < 0.000001)
                    {
                        if (b.ux)
                            BCsIndex.Add(6 * id);

                        if (b.uy)
                            BCsIndex.Add(6 * id + 1);

                        if (b.uz)
                            BCsIndex.Add(6 * id + 2);

                        if (b.rx)
                            BCsIndex.Add(6 * id + 3);

                        if (b.ry)
                            BCsIndex.Add(6 * id + 4);

                        if (b.rz)
                            BCsIndex.Add(6 * id + 5);
                    }
                }

            }

            return BCsIndex;
        }


        private Vector<double> CreateLoadList(List<LoadClass> _lc, List<NodeClass> _nodes)
        {

            Vector<double> LoadValue = DenseVector.OfEnumerable(new double[_nodes.Count * 6]);

            foreach (var load in _lc)

            {
                List<Point3d> LoadPts = new List<Point3d>();

                if (load.Id == true)
                {
                    LoadPts.Add(load.coordinate);
                }


                else
                {
                    LoadPts.Add(load.stPt);
                    LoadPts.Add(load.enPt);
                }


                Vector3d LoadVec = load.LoadVec;

                for (int i = 0; i < _nodes.Count; i++)
                {
                    Point3d node = _nodes[i].pt;
                    int id = _nodes[i].Id;


                    for (int j = 0; j < LoadPts.Count; j++)
                    {
                        Point3d loadLocation = LoadPts[j];
                        if (loadLocation.DistanceTo(node) < 0.00001)
                        {
                            if (load.Id1)
                            {
                                LoadValue[id * 6] += LoadVec.X;
                                LoadValue[id * 6 + 1] += LoadVec.Y;
                                LoadValue[id * 6 + 2] += LoadVec.Z;
                            }

                            else
                            {
                                LoadValue[id * 6 + 3] += load.mVal[0];
                                LoadValue[id * 6 + 4] += load.mVal[1];
                                LoadValue[id * 6 + 5] += load.mVal[2];
                            }

                        }
                        else
                        {
                            LoadValue[id * 6] += 0;
                            LoadValue[id * 6 + 1] += 0;
                            LoadValue[id * 6 + 2] += 0;
                            LoadValue[id * 6 + 3] += 0;
                            LoadValue[id * 6 + 4] += 0;
                            LoadValue[id * 6 + 5] += 0;
                        }

                    }

                }

            }

            return LoadValue;

        }

        private static void CreateForces(List<BeamClassLinear> bars, List<Point3d> points, Vector<double> _def, List<Matrix<double>> Lk_eg, Matrix<double> N, Matrix<double> dN, List<Matrix<double>> T_List, List<Matrix<double>> T_List_six, out Vector<double> forces, out Vector<double> moment, out Vector<double> strain, out Vector<double> stress, out List<Vector<double>> _u, out double M_max, out double umax, out List<Line> crv)
        {
            //Matrix<double> k_eG = DenseMatrix.OfArray(new double[6, 6]);
            Vector<double> u = DenseVector.OfEnumerable(new double[6]);
            Vector<double> eps = DenseVector.OfEnumerable(new double[6]);
            Vector<double> sigma = DenseVector.OfEnumerable(new double[6]);
            Vector<double> S;

            Vector<double> v = DenseVector.OfEnumerable(new double[12]);
            Vector<double> vl = DenseVector.OfEnumerable(new double[12]);
            var _M = new double();

            Vector<double> mom = DenseVector.OfEnumerable(new double[points.Count * 3]);
            Vector<double> forc = DenseVector.OfEnumerable(new double[points.Count * 3]);
            Vector<double> dNN = DenseVector.OfEnumerable(new double[12]);
            
            List<Line> bending = new List<Line>();

            var u_max = new double();
            double u_Max = 0;

            var m_max = new double();
            double m_Max = 0;


            List<Vector<double>> u_lst = new List<Vector<double>>();
            List<double> M_lst = new List<double>();

            foreach (BeamClassLinear b in bars)
            {


                var E = b.material.youngsModolus;
                var my = b.material.v;
                int node1 = b.startNode.Id;
                int node2 = b.endNode.Id;

                v[0] = _def[node1 * 6];
                v[1] = _def[node1 * 6 + 1];
                v[2] = _def[node1 * 6 + 2];
                v[3] = _def[node1 * 6 + 3];
                v[4] = _def[node1 * 6 + 4];
                v[5] = _def[node1 * 6 + 5];
                v[6] = _def[node2 * 6];
                v[7] = _def[node2 * 6 + 1];
                v[8] = _def[node2 * 6 + 2];
                v[9] = _def[node2 * 6 + 3];
                v[10] = _def[node2 * 6 + 4];
                v[11] = _def[node2 * 6 + 5];


                //displacement

                var B = dN;
                Matrix<double> T_six_trans = T_List_six[b.Id].Transpose();
                vl = T_List[b.Id] * v;       //transponerer v

                u = N.Multiply(vl);
                u = T_six_trans.Multiply(u);


                u_lst.Add(u);

                u_max = Math.Sqrt(Math.Pow(u[0], 2) + Math.Pow(u[1], 2) + Math.Pow(u[2], 2));

                if (u_max > u_Max)
                {
                    u_Max = u_max;
                }

                Curve currentLine = b.axis;
                double L = Math.Round(currentLine.GetLength() * 1000.00, 2);
                List<Point3d> curve_pts = new List<Point3d>();
                double n = 2;
                var x = 0.0;

                for (int i = 0; i < n + 1; i++)

                {
                    double ddN3 = -6.0 / Math.Pow(L, 2) + 12.0 * x / Math.Pow(L, 3);
                    double ddN4 = -6.0 * x / Math.Pow(L, 2) + 4.0 / L;
                    double ddN5 = 6.0 / Math.Pow(L, 2) - 12.0 * x / Math.Pow(L, 3);
                    double ddN6 = -6.0 * x / Math.Pow(L, 2) + 2.0 / L;

                    dNN = DenseVector.OfEnumerable(new double[]

                    {0,     0,   ddN3,  0,   ddN4,      0,       0,     0,    ddN5,    0,   ddN6,    0}
                    );

                    double wxx = dNN.DotProduct(vl);
                    _M = E * b.section.Iy * wxx;       //langs xy plan. varierer langs z ( gange med h/2)
                    M_lst.Add(_M);

                    x += L / n;

                    Vector3d vector = b.axis.CurvatureAt((L / n) / L);
                    Vector3d nullVec = new Vector3d(0,0,0);
                    vector.Unitize();
                    Point3d st = b.axis.PointAtNormalizedLength((L / n * i) / L);
                    
                    if (vector == nullVec)
                    {
                        vector = new Vector3d(0,0,1);
                    }
                    
                    Point3d en = new Point3d(st.X - vector.X * _M / 10000000.0, st.Y - vector.Y * _M / 10000000.0, st.Z - vector.Z * _M / 10000000.0);

                    
                    curve_pts.Add(en);

                    
                }



                eps = B.Multiply(vl);        //strain (3)



                sigma = E * eps;            //stress = E*3

                S = Lk_eg[b.Id].Multiply(v);

                forc[node1 * 3] += S[0];
                forc[node1 * 3 + 1] += S[1];
                forc[node1 * 3 + 2] += S[2];

                forc[node2 * 3] += S[6];
                forc[node2 * 3 + 1] += S[7];
                forc[node2 * 3 + 2] += S[8];


                mom[node1 * 3] += S[3];
                mom[node1 * 3 + 1] += S[4];
                mom[node1 * 3 + 2] += S[5];

                mom[node2 * 3] += S[9];
                mom[node2 * 3 + 1] += S[10];
                mom[node2 * 3 + 2] += S[11];



                for (int a = 1; a < curve_pts.Count; a++)
                {
                    Line line = new Line(curve_pts[a - 1], curve_pts[a]);
                    bending.Add(line);
                }
                Line line1 = new Line(b.startNode.pt, curve_pts[0]);
                Line line2 = new Line(b.endNode.pt, curve_pts[curve_pts.Count - 1]);
                bending.Add(line1);
                bending.Add(line2);


            }


            


            foreach (double m_temp in M_lst)
            {
                
                if (Math.Abs(m_temp) > Math.Abs(m_Max))
                {
                    m_Max = m_temp;
                }
            }



            forces = forc;
            moment = mom;
            strain = eps;
            stress = sigma;
            umax = u_Max;
            _u = u_lst;
            M_max = m_Max;
            crv = bending;

        }


        private static void CreateGenerelizedShapeFunc(List<BeamClassLinear> bars, Matrix<double> k_tot, double X, out Matrix<double> N, out Matrix<double> dN)
        {

            Matrix<double> _N = DenseMatrix.OfArray(new double[6, 12]);
            Matrix<double> _dN = DenseMatrix.OfArray(new double[6, 12]);



            foreach (BeamClassLinear b in bars)
            {
                //trying






                Curve currentLine = b.axis;
                double L = currentLine.GetLength();
                double x = X * L;


                double N1 = 1.0 - x / L;                                                                      //axial translation in node 1
                double N2 = x / L;                                                                          //axial translation in node 2
                double N3 = 1.0 - 3.00 * Math.Pow(x, 2) / Math.Pow(L, 2) + 2.0 * Math.Pow(x, 3) / Math.Pow(L, 3);  //translation node 1
                double N4 = -x * (1.0 - 2.00 * x / L + Math.Pow(x, 2) / Math.Pow(L, 2));                                 //rotation node 1
                double N5 = 3.00 * Math.Pow(x, 2) / Math.Pow(L, 2) - 2.0 * Math.Pow(x, 3) / Math.Pow(L, 3);      //translation node 2
                double N6 = x * (x / L - Math.Pow(x, 2) / Math.Pow(L, 2));                                 //rotation node 2


                double dN1 = -1.0 / L;
                double dN2 = 1.0 / L;
                double dN3 = -6.0 * x / Math.Pow(L, 2) + 6.0 * Math.Pow(x, 2) / Math.Pow(L, 3);
                double dN4 = -3.0 * Math.Pow(x, 2) / Math.Pow(L, 2) + 4 * x / L - 1.0;
                double dN5 = -dN3;
                double dN6 = -3.0 * Math.Pow(x, 2) / Math.Pow(L, 2) + 2.0 * x / L;

                double ddN3 = -6.0 / Math.Pow(L, 2) + 12.0 * x / Math.Pow(L, 3);
                double ddN4 = -6.0 * x / Math.Pow(L, 2) + 4.0 / L;
                double ddN5 = 6.0 / Math.Pow(L, 2) - 12.0 * x / Math.Pow(L, 3);
                double ddN6 = -6.0 * x / Math.Pow(L, 2) + 2.0 / L;


                // første rad blir ux = N1*ux1 + N2*ux2
                //andre rad blir uy = N3*uy1 + N5*uy2 etc.

                _N = Matrix<double>.Build.DenseOfArray(new double[,]
                {
                {N1,   0,     0,    0,    0,     0,     N2,     0,     0,     0,     0,    0},
                {0,   N3,     0,    0,    0,    -N4,     0,    N5,     0,     0,     0,   -N6},
                {0,    0,    N3,    0,   N4,     0,      0,     0,    N5,     0,   N6,    0},
                {0,    0,     0,    N1,   0,     0,      0,     0,     0,    N2,     0,    0},
                {0,    0,    dN3,   0,  dN4,     0,      0,     0,    dN5,    0,   dN6,   0},
                {0,   dN3,    0,    0,    0,    -dN4,    0,    dN5,    0,     0,     0,   -dN6},
                

                });



                //N1 and N2 is linear shape function and is derived one (Translationin x dir) while the rest is derived twice for beams. 

                _dN = Matrix<double>.Build.DenseOfArray(new double[,]
                {
                {dN1,   0,    0,    0,     0,       0,    dN2,     0,      0,     0,     0,      0},
                {0,    dN3,   0,    0,     0,      -dN4,     0,    dN5,     0,     0,     0,    -dN6},
                {0,     0,   dN3,   0,   dN4,       0,       0,     0,    dN5,     0,   dN6,     0},
                {0,     0,    0,   dN1,    0,       0,       0,     0,      0,    dN2,    0,      0},
                {0,     0,   ddN3,  0,   ddN4,      0,       0,     0,    ddN5,    0,   ddN6,    0},
                {0,    ddN3,  0,    0,     0,     -ddN4,     0,    ddN5,    0,     0,     0,   -ddN6},
                

                });



            }

            N = _N;
            dN = _dN;
        }






        private static void CreateGlobalStiffnesMatrix(List<BeamClassLinear> bars, List<Point3d> points, out Matrix<double> k_tot, out List<Matrix<double>> L_k_eg, out List<Matrix<double>> T_list, out List<Matrix<double>> T_list_six)

        {

            int dofs = points.Count * 6;
            Matrix<double> K_tot = DenseMatrix.OfArray(new double[dofs, dofs]);
            Matrix<double> K_eG = DenseMatrix.OfArray(new double[12, 12]);
            List<Matrix<double>> LK_eG = new List<Matrix<double>>();

            List<Matrix<double>> transformationmatrix = new List<Matrix<double>>();
            List<Matrix<double>> tranmatrix_six = new List<Matrix<double>>();

            foreach (BeamClassLinear b in bars)
            {


                Curve currentLine = b.axis;
                double L = currentLine.GetLength() * 1000;



                double X = b.section.CSA * b.material.youngsModolus / L;
                double Y1 = 12.00 * b.material.youngsModolus * b.section.Iz / (Math.Pow(L, 3));
                double Y2 = 6.00 * b.material.youngsModolus * b.section.Iz / (Math.Pow(L, 2));  //slender profile, approx.
                double Y3 = 4.00 * b.material.youngsModolus * b.section.Iz / L;
                double Y4 = 2.00 * b.material.youngsModolus * b.section.Iz / L;

                double Z1 = 12.00 * b.material.youngsModolus * b.section.Iy / (Math.Pow(L, 3));
                double Z2 = 6.00 * b.material.youngsModolus * b.section.Iy / (Math.Pow(L, 2));  //slender profile, approx.
                double Z3 = 4.00 * b.material.youngsModolus * b.section.Iy / L;
                double Z4 = 2.00 * b.material.youngsModolus * b.section.Iy / L;

                double C = (b.material.G * b.section.J) / L;


                Point3d p1 = new Point3d(Math.Round(currentLine.PointAtStart.X, 6), Math.Round(currentLine.PointAtStart.Y, 6), Math.Round(currentLine.PointAtStart.Z, 6));
                Point3d p2 = new Point3d(Math.Round(currentLine.PointAtEnd.X, 6), Math.Round(currentLine.PointAtEnd.Y, 6), Math.Round(currentLine.PointAtEnd.Z, 6)); 

                Vector3d unitX = new Vector3d(1, 0, 0);
                Vector3d unitY = new Vector3d(0, 1, 0);
                Vector3d unitZ = new Vector3d(0, 0, 1);


                Vector3d nullvec = new Vector3d(0, 0, 0);
                Vector3d vec1x = b.axis.TangentAt(0);
                Vector3d vec1y = new Vector3d(-vec1x.Y, vec1x.X, 0);
                Vector3d vec1z = Vector3d.CrossProduct(vec1y, vec1x);

                vec1z.Unitize();
                vec1x.Unitize();
                vec1y.Unitize();
                
                double _1l1 = Math.Cos(Vector3d.VectorAngle(vec1x, unitX));
                double _1m1 = Math.Cos(Vector3d.VectorAngle(vec1x, unitY));
                double _1n1 = Math.Cos(Vector3d.VectorAngle(vec1x, unitZ));

                double _1l2 = Math.Cos(Vector3d.VectorAngle(vec1y, unitX));
                double _1m2 = Math.Cos(Vector3d.VectorAngle(vec1y, unitY));
                double _1n2 = Math.Cos(Vector3d.VectorAngle(vec1y, unitZ));

                double _1l3 = Math.Cos(Vector3d.VectorAngle(vec1z, unitX));
                double _1m3 = Math.Cos(Vector3d.VectorAngle(vec1z, unitY));
                double _1n3 = Math.Cos(Vector3d.VectorAngle(vec1z, unitZ));

                Matrix<double> t = DenseMatrix.OfArray(new double[,]
                {
                        {_1l1,   _1m1,     _1n1},
                        {_1l2,   _1m2,     _1n2},
                        {_1l3,   _1m3,     _1n3},
                });


                t.CoerceZero(0.00001);

                var T_t = t.DiagonalStack(t);
                var T = T_t.DiagonalStack(T_t);

                transformationmatrix.Add(T);
                tranmatrix_six.Add(T_t);



                Matrix<double> ke = DenseMatrix.OfArray(new double[,]
                {
                        {X,  0,   0,   0,   0,   0,  -X,   0,   0,   0,   0,   0},
                        {0,  Y1,  0,   0,   0,   Y2,  0,  -Y1,  0,   0,   0,   Y2},
                        {0,  0,   Z1,  0,  -Z2,  0,   0,   0,  -Z1,  0,  -Z2,  0},
                        {0,  0,   0,   C,   0,   0,   0,   0,   0,  -C,   0,   0},
                        {0,  0 , -Z2,  0,   Z3,  0,   0,   0,   Z2,  0,   Z4,  0},
                        {0,  Y2,  0,   0,   0,   Y3,  0,  -Y2,  0,   0,   0,   Y4},
                        {-X, 0,   0,   0,   0,   0,   X,   0,   0,   0,   0,   0},
                        {0, -Y1,  0,   0,   0,  -Y2,  0,   Y1,  0,   0,   0,  -Y2},
                        {0,  0,  -Z1,  0,   Z2,  0,   0,   0,   Z1,  0,   Z2,  0},
                        {0,  0,   0,  -C,   0,   0,   0,   0,   0,   C,   0,   0},
                        {0,  0,  -Z2,  0,   Z4,  0,   0,   0,   Z2,  0,   Z3,  0},
                        {0,  Y2,  0,   0,   0,   Y4,  0,  -Y2,  0,   0,   0,   Y3},


                });
                Matrix<double> Tt = T.Transpose(); //transpose
                Matrix<double> KG = Tt.Multiply(ke);
                K_eG = KG.Multiply(T);

                K_eG.CoerceZero(0.00001);

                LK_eG.Add(K_eG);

                int node1 = b.startNode.Id;
                int node2 = b.endNode.Id;

                for (int i = 0; i < K_eG.RowCount / 2; i++)
                {
                    for (int j = 0; j < K_eG.ColumnCount / 2; j++)
                    {
                        K_tot[node1 * 6 + i, node1 * 6 + j] += Math.Round(K_eG[i, j], 5);
                        K_tot[node1 * 6 + i, node2 * 6 + j] += Math.Round(K_eG[i, j + 6], 5);
                        K_tot[node2 * 6 + i, node1 * 6 + j] += Math.Round(K_eG[i + 6, j], 5);
                        K_tot[node2 * 6 + i, node2 * 6 + j] += Math.Round(K_eG[i + 6, j + 6], 5);

                    }

                }


            }

            k_tot = K_tot;
            L_k_eg = LK_eG;
            T_list = transformationmatrix;
            T_list_six = tranmatrix_six;


        }

        private static void CreateReducedStiffnesMatrix(List<int> _BC, Matrix<double> K_tott, Vector<double> Rvec, out Matrix<double> K_red)
        {


            for (int i = 0; i < _BC.Count; i++)
            {
                int sizeOfMatrix = K_tott.ColumnCount;
                int blockedIndex = _BC[i];

                Rvec[blockedIndex] = 0;

                for (int r = 0; r < sizeOfMatrix; r++)
                {
                    K_tott[blockedIndex, r] = 0;
                    K_tott[r, blockedIndex] = 0;

                    if (K_tott[r, r] == 0)
                    {
                        K_tott[r, r] = 1;
                    }
                }
            }

            K_red = K_tott;



        }


        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                return Master.Properties.Resources.Linear.ToBitmap();
                //return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("e989e816-9e07-4039-bff9-209a2b3f8a2c"); }
        }
    }
}