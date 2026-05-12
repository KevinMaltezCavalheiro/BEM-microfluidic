int compute_G_matrix(map<string, pair<int, int>> &groups,vector<size_t> &all_nodeTags,vector<double> &all_nodecoord,vector<double> &all_nodeparametricCoord,vector<string> &keys, string &groupname,vector<int> &tags,string &integration_rule, vectorpair &dimTags)
{
    //le but est d'avoir un vecteur qui contient tous les éléments de la frontière. Pour ce faire, on combine tous les physical groups (right_edge+top_edge+left_edge+bottom_edge dans le cas du rectangle)
    vector<double> boundary_elements_tags(0,0);
    vector<double> boundary_sort_nodes_tags(0,0);
    vector<double> boundary_nodes_coordinates_x(0,0);
    vector<double> boundary_nodes_coordinates_y(0,0);

    vector<double> coord;
    vector<double> parametricCoord;
    vector<size_t> nodeTags;

    model::mesh::renumberNodes();
    model::mesh::getNodes(nodeTags, coord, parametricCoord);

    // display the vectors
    cout << "\tnodeTags = ";
    for(int i=0; i<nodeTags.size(); ++i)
    {
        cout << nodeTags[i] << " ";
        boundary_nodes_coordinates_x.push_back(coord[i*3]);//coordonnée x de tous les noeuds de la frontière
        boundary_nodes_coordinates_y.push_back(coord[(i*3)+1]);//coordonnée y de tous les noeuds de la frontière
    }
    cout << endl;

    //just for check
    for(int i=0; i<nodeTags.size(); ++i)
    {
        cout << nodeTags[i] << "(" << coord[i*3] << "," << coord[(i*3)+1] << ") ";
    }
    cout << endl;

    cout << "\tcoord = ";
    for(int i=0; i<coord.size(); ++i)
    {
        cout << coord[i] << " ";
    }
    cout << endl;

    //    //le problème majeur, c'est qu'on a un gros vecteur, avec tous les noeuds sauf que leur numérotation n'est pas la même que le sens de parcours de la boundary (on aimerait bien qu'ils se suivent). (typiquement aller deflag fltk::run(); aller dans tools, option, mesh, nodes label)
    //    boundary_sort_nodes_tags.resize(nodeTags.size(), 0);
    //    vector<double> boundary_sort_length_elements(nodeTags.size(),0);
    //    sorting_nodes(nodeTags,boundary_nodes_coordinates_x,boundary_nodes_coordinates_y,boundary_sort_nodes_tags,boundary_sort_length_elements);

    //code typique de parcours des éléments: première boucles, les différents physicals groups (left_edge,right_edge,...)
    for (int i = 0; i < dimTags.size(); ++i)
    {
        int dim = dimTags[i].first;
        int tag = dimTags[i].second;
        model::getEntitiesForPhysicalGroup(dim, tag, tags);
        string name;
        model::getPhysicalName(dim, tag, name);
        groups[name] = {dim, tag};
        cout << endl;
        cout << endl;
        cout << "Entities in group named: " << name << endl;

        //deuxième boucle: loop over entities (au sein d'un physical group) par exemple j'aurais pu créer un seul physical group boundary et au sein de se groupe, cette boucle aurait parcouru les 4 lignes qui forment les différents bords (que j'aurais regroupé dans boundary formée de 4 lignes) (ici on n'a que une courbe (un bord) par groupe, on aurait pu en mettre plusieurs)
        for(int k = 0; k < tags.size(); ++k)
        {
            cout << "Entity (" << dim << "," << tags[k] << "):\n";
            vector<int> elementTypes;
            vector<vector<size_t>> elementTags;
            vector<vector<size_t>> nodeTagsk;

            model::mesh::getElements(elementTypes, elementTags, nodeTagsk, dim, tags[k]);

            // loop over element types ([elements] and [nodes]) looping over the elements means looping over each element type of the mesh of each entity of a given physical group. Each entity can have several sets of elements of different types (in our example, we only have 1 set of triangles). par exemple an aurait pu avoir des carrés puis des triangles
            for (int l = 0; l < elementTypes.size(); ++l)
            {
                //loop sur un type d'élement par exemple sur tous les triangles
                //just to check
                cout << "\telement type " << elementTypes[l] << ":\n";
                cout << "\t\telements: ";
                for (int j = 0; j < elementTags[l].size(); ++j)
                {
                    cout << elementTags[l][j] << " ";
                }
                cout << endl;
                cout << "\t\tnodes: ";
                for (int j = 0; j < nodeTagsk[l].size(); ++j)
                {
                    cout << nodeTagsk[l][j] << "(" << coord[nodeTagsk[l][j]] << ") ";
                }
            }
        }
    }

    cout << endl;
    return 0;
}

void sorting_nodes(vector<size_t> &nodeTags,vector<double> &boundary_nodes_coordinates_x,vector<double> &boundary_nodes_coordinates_y,vector<double> &boundary_sort_nodes_tags, vector<double> &boundary_sort_length_elements)
{
    double length = 0; //longueure entre 2 noeuds à minimser
    double temp = 0; //variable temporaire pour éviter de stupides erreurs
    double min_tag_prev = 0; //retient le tag du précédent noeud, on pourrait s'en passer mais c'est plus simple pour raisonner
    double min_tag_next = 0;//retient le tag du noeud suivant
    double min_value = 0;//valeur minimale de la longueur pour repérer 2 noeuds successifs
    boundary_sort_nodes_tags[0] = 1;//on initialise le premier élément (là où l'on commence la boucle des noeuds puisque la frontière est fermée) ici à l'élément 1.
    bool already_count = 0;//si le noeud minimum fait partie des noeuds déjà considérés, on ne le compte pas et on passe au suivant

    //parcours le vecteur des tags triés pour le complèter. Les deux premiers éléments ont déjà été calculé.
    for(int i = 1; i < nodeTags.size(); i++)
    {
        min_value = HUGE_VAL; //pour chercher le minimum, on démarre d'une valeur très grande
        //on cherche parmis tous les autres noeuds lequel est le minimum. Ils ne sont pas dans l'ordre, on doit tous les considérer.
        for(int j = 0; j < nodeTags.size(); j++)
        {
            already_count = 0;// initialement le noeud n'est pas compté
            //on parcours tous les éléments du vecteur node_tag_sort pour vérifier que le noeud potentiellement minimum j n'en fait pas partie. Si c'est le cas, break et on passe au noeud suivant
            for(int k = 0; k < i; k++)
            {
                if(j+1 == boundary_sort_nodes_tags[k])
                {
                    already_count = 1;
                    break;
                }
            }

            // si already count == 1, le noeud est déjà pris en compte et on passe au j suivant
            if(already_count == 0)
            {
                //variable temporaire: (x_next-x_previous)^2+(y_next-y_previous)^2
                temp = (boundary_nodes_coordinates_x[j]-boundary_nodes_coordinates_x[min_tag_prev])*(boundary_nodes_coordinates_x[j]-boundary_nodes_coordinates_x[min_tag_prev]) + (boundary_nodes_coordinates_y[j]-boundary_nodes_coordinates_y[min_tag_prev])*(boundary_nodes_coordinates_y[j]-boundary_nodes_coordinates_y[min_tag_prev]);
                length = sqrt(temp); //calcule la longueure

                if(length <= min_value)
                {
                    min_tag_next = j;
                    min_value = length;
                }
            }
        }
        //just for check
        //cout << "min tag prev: " << min_tag_prev+1 << endl;
        boundary_sort_nodes_tags[i] = min_tag_next+1;
        boundary_sort_length_elements[i] = min_value;
        min_tag_prev = min_tag_next;
        //        cout << "min tag next: " << min_tag_next+1 << endl;
        //        cout << "min value: " << min_value << endl;
        //        cout << endl;
    }
    length = 0;
    int first = boundary_sort_nodes_tags.front()-1;
    int last = boundary_sort_nodes_tags.back()-1;
    temp = (boundary_nodes_coordinates_x[first]-boundary_nodes_coordinates_x[last])*(boundary_nodes_coordinates_x[first]-boundary_nodes_coordinates_x[last]) + (boundary_nodes_coordinates_y[first]-boundary_nodes_coordinates_y[last])*(boundary_nodes_coordinates_y[first]-boundary_nodes_coordinates_y[last]);
    length = sqrt(temp); //calcule la longueure
    boundary_sort_length_elements[0] = length;//pour terminer la boucle, on relie le premier au dernier élément
}

//void get_length_elements(vector<double> &boundary_sort_nodes_tags, vector<double> &boundary_length_elements, vector<double> &boundary_nodes_coordinates_x,vector<double> &boundary_nodes_coordinates_y)
//{
//    for(int i = 0; i < boundary_sort_nodes_tags.size()-1; i++)
//    {
//        boundary_length_elements[i] = sqrt((boundary_nodes_coordinates_x[boundary_sort_nodes_tags[i+1]]-boundary_nodes_coordinates_x[boundary_sort_nodes_tags[i]])*(boundary_nodes_coordinates_x[boundary_sort_nodes_tags[i+1]]-boundary_nodes_coordinates_x[boundary_sort_nodes_tags[i]])+(boundary_nodes_coordinates_y[boundary_sort_nodes_tags[i+1]]-boundary_nodes_coordinates_y[boundary_sort_nodes_tags[i]])*(boundary_nodes_coordinates_y[boundary_sort_nodes_tags[i+1]]-boundary_nodes_coordinates_y[boundary_sort_nodes_tags[i]]));
//    }
//    boundary_length_elements[boundary_length_elements.size()-1] = sqrt((boundary_nodes_coordinates_x[boundary_sort_nodes_tags[boundary_length_elements.size()-1]]-boundary_nodes_coordinates_x[boundary_sort_nodes_tags[0]])*(boundary_nodes_coordinates_x[boundary_sort_nodes_tags[boundary_length_elements.size()-1]]-boundary_nodes_coordinates_x[boundary_sort_nodes_tags[0]])+(boundary_nodes_coordinates_y[boundary_sort_nodes_tags[boundary_length_elements.size()-1]]-boundary_nodes_coordinates_y[boundary_sort_nodes_tags[0]])*(boundary_nodes_coordinates_y[boundary_sort_nodes_tags[boundary_length_elements.size()-1]]-boundary_nodes_coordinates_y[boundary_sort_nodes_tags[0]]));
//}


//int boundary_conditions()
//{
//    /*-----------get physical properties of boundary----------------------*/
//    onelab::getNames(keys, "(Boundary Conditions).+");
//    //exemple: "Boundary Conditions/bottom_edge/tx"
//    for (auto &key : keys)
//    {
//        // get corresponding value
//        vector<double> value;
//        onelab::getNumber(key, value);
//        // expected key structure is "type/group_name/field"
//        // => split the key string into components
//        stringstream ss(key);
//        vector<string> words;
//        string word;
//        while(getline(ss, word, '/'))// read string until '/'
//        {
//            words.push_back(word);
//        }
//        //words[0] = Boundary Conditions, words[1] = bottom_edge, words[2] = tx
//        if(words.size()==3)
//        {
//            groupname = ""; // will remain empty only if tx == ty == 0
//            if(words[2].compare("tx") == 0)
//            {
//                tx = value[0];
//                string new_key = "Boundary Conditions/" + words[1] + "/ty";
//                onelab::getNumber(new_key, value);
//                ty = value[0];
//                if(tx != 0 || ty != 0) // there is no need the consider Neumann BC if it is homogeneous
//                {
//                    groupname = words[1];
//                }
//            }
//            //here, if we are dealing with neumann BC's, we should have retrieved a non-empty group_name
//            if(groupname.compare("") != 0)
//            {
//                cout << "------------------------------------" << "\n";
//                cout << groupname << ", tx: " << tx << ", ty: " << ty << "\n";
//                //loop over elements of group_name
//                int dim = groups[groupname].first;
//                int tag = groups[groupname].second;
//                if (tag == 0)
//                {
//                    cerr << "Group '" << groupname << "' does not exist!\n";
//                    return 1;
//                }
//                model::getEntitiesForPhysicalGroup(dim, tag, tags);
//
//                cout << "Entities in group named '" << groupname << "': ";
//                for (int i = 0; i < tags.size(); ++i)
//                {
//                    cout << tags[i] << " ";
//                }
//                cout << endl;
//                // loop over entities
//                for (int k = 0; k < tags.size(); ++k)
//                {
//                    cout << "Entity (" << dim << "," << tags[k] << "):" << endl;
//                    vector<int> elementTypes;
//                    vector<vector<size_t>> elementTags;
//                    vector<vector<size_t>> nodeTags;
//                    model::mesh::getElements(elementTypes, elementTags, nodeTags, dim, tags[k]);
//                    //looping over element types
//                    for (int i = 0; i < elementTypes.size(); ++i)
//                    {
//                        //just to check
//                        cout << "\telement type " << elementTypes[i] << ":" << endl;
//                        cout << "\t\telements: ";
//                        for (int j = 0; j < elementTags[i].size(); ++j)
//                        {
//                            cout << elementTags[i][j] << " ";
//                        }
//                        cout << endl;
//                        cout << "\t\tnodes: ";
//                        for (int j = 0; j < nodeTags[i].size(); ++j)
//                        {
//                            cout << nodeTags[i][j] << " ";
//                        }
//                        cout << endl;
//
//                        // getting element properties
//                        string elementName;
//                        int dim, order, numNodes, numPrimaryNodes;
//                        vector<double> localNodeCoord;
//                        model::mesh::getElementProperties(elementTypes[i],
//                                                          elementName, dim, order,
//                                                          numNodes, localNodeCoord,
//                                                          numPrimaryNodes);
//
//                        // get Gauss points coordinates and weights
//                        vector<double> localCoords, weights;
//                        model::mesh::getIntegrationPoints(elementTypes[i], integration_rule, localCoords, weights);//maybe use different integration method on boundary than in domain
//                        cout << "\t" << weights.size() << " integration points" << endl;
//
//                        cout << "\tweights: ";
//                        for (int j = 0; j < weights.size(); ++j)
//                        {
//                            cout << weights[j] << " ";
//                        }
//                        cout << endl;
//
//                        cout << "\tlocalCoords (gmsh): ";
//                        for (int j = 0; j < localCoords.size(); ++j)
//                        {
//                            cout << localCoords[j] << " ";
//                        }
//                        cout << endl;
//
//                        // get determinants and Jacobians (only determinants will be useful)
//                        vector<double> jacobians, determinants, coords;
//                        model::mesh::getJacobians(elementTypes[i], localCoords, jacobians, determinants, coords,tags[k]); //ATTENTION ai dû rajouter mon tags[k] pr focus on this particular entity
//
//                        //just to check
//                        cout << "\tGot " << determinants.size()
//                        << " Jacobians for type " << elementTypes[i] << endl;
//
//                        cout << "\tdeterminants (gmsh): ";
//                        for (int j = 0; j < determinants.size(); ++j)
//                        {
//                            cout << determinants[j] << " ";
//                        }
//                        cout << endl;
//
//                        int numComponents, numOrientations;
//                        vector<double> basisFunctions;
//                        vector<double> basisFunctionsGradient;
//                        model::mesh::getBasisFunctions(elementTypes[i], localCoords, "Lagrange", numComponents, basisFunctions, numOrientations);
//
//                        //just to check
//                        cout << "\tGot " << basisFunctions.size() << " basis functions for type " << elementTypes[i] << endl;
//                        for (int j = 0; j < basisFunctions.size(); ++j)
//                        {
//                            cout << basisFunctions[j] << " ";
//                        }
//                        cout << endl;
//
//                        model::mesh::getBasisFunctions(elementTypes[i], localCoords,
//                                                       "GradLagrange", numComponents,
//                                                       basisFunctionsGradient,
//                                                       numOrientations);
//
//                        /*-----------computing the full_G_matrix---------*/
//
//                        for (int el = 0; el < elementTags[i].size(); el++)
//                        {
//                            //for (int el = 0; el < 1; el++){
//                            // elemental K matrix and f vector initialization
//                            Matrix<double, Dynamic, Dynamic> G_matrix(numNodes,numNodes);
//                            vector<double> f_vector(numNodes);
//                            for (int j = 0; j < numNodes; ++j)
//                            {
//                                for (int l = 0; l < numNodes; ++l)
//                                {
//                                    G_matrix(j,l)=0.;
//                                }
//                            }
//
//                            for(int j = 0; j < weights.size(); ++j)
//                            {
//                                // converting jacobian to Eigen format
//                                Map<Matrix<double, Dynamic, Dynamic> >
//                                jacob(&jacobians[9*j+9*weights.size()*el], 3, 3);
//                                // 3D to 2D (plane stress case)
//                                Matrix<double, 2, 2> jacob2D;
//                                jacob2D(0,0) = jacob(0,0);
//                                jacob2D(1,0) = jacob(1,0);
//                                jacob2D(0,1) = jacob(0,1);
//                                jacob2D(1,1) = jacob(1,1);
//
//                                // computing the inverse
//                                Matrix<double, 2, 2> jacobinv = jacob2D.inverse();
//
//                                // computing the transpose of the inverse
//                                Matrix<double, 2, 2> jacobinvtrans = jacobinv.transpose();
//                                //compute the G_Matrix
//
//                                for(int l = 0; l < numNodes; ++l)
//                                {// looping over the number of shape functions
//                                    for(int m = 0; m < numNodes; ++m)
//                                    {
//                                        if(m != l)
//                                        {
//
//                                        }
//                                        if(m == l)
//                                        {
//
//                                        }
//                                    }
//                                }
//                            }
//                            /*-----------computing the f vector components related to this Neumann BC---------*/
//                            /*//looping over the elements
//                             vector<double> f_vector(2*numNodes);
//                             for (int el = 0; el < elementTags[i].size(); el++){
//                             for (int j = 0; j < 2*numNodes; ++j){
//                             f_vector[j] = 0.;
//                             }
//                             // looping over the Gauss points (formula of slide 17)
//                             for (int j = 0; j < weights.size(); ++j){
//                             // looping over the number of shape functions
//                             for (int l = 0; l < numNodes; ++l){
//                             f_vector[2*l] = f_vector[2*l] + basisFunctions[l+numNodes*j]*tx*determinants[j+weights.size()*el]*weights[j]; //pas oublier le dét !
//                             f_vector[2*l+1] = f_vector[2*l+1] + basisFunctions[l+numNodes*j]*ty*determinants[j+weights.size()*el]*weights[j]; // le weight non plus !
//                             }
//                             }
//                             cout << "\telemental f vector of element " << elementTags[i][el] << ": ";
//                             for (int j = 0; j < f_vector.size(); ++j)
//                             cout << f_vector[j] << " ";
//                             cout << '\n';*/
//                            /*------------- ASSEMBLY PROCESS ------------*/
//                            /*for (int j = 0; j < 2*numNodes; ++j)
//                             {
//                             int first_node_index = (int) j/2;
//                             int first_node = nodeTags[i][numNodes*el+first_node_index];
//                             int my_row = 2*(first_node-1) + j % 2; // node list begins at 1, '-1' correction to begin at 0
//                             full_f_vector[my_row] += f_vector[j];
//                             //cout << my_row << "\n";
//                             }
//                             }*/
//                        }
//                    }
//                }
//            }
//        }
//    }
//    return 0;
//}


//int initialise_G_matrix(vector<size_t> &all_nodeTags, map<string, pair<int, int>> &groups, string &integration_rule, double bx, double rho, double by, Matrix<double, 3, 3> &G_matrix, Matrix<double, 3, 3> &full_G_matrix, string groupname)
//{
//    //attribue dim et tag au groupname = "boundary". attention, il doit map avec le fichier geo
//    int dim = groups[groupname].first;
//    int tag = groups[groupname].second;
//    //on s'assure que le vecteur n'est pas vide
//    if (tag == 0)
//    {
//        cerr << "Group '" << groupname << "' does not exist!\n";
//        return 1;
//    }
//
//    // getting entities of the boundary
//    vector<int> tags;
//    model::getEntitiesForPhysicalGroup(dim, tag, tags);
//
//    cout << "Entities in group named '" << groupname << "': ";
//    for (int i = 0; i < tags.size(); ++i)
//    {
//        cout << tags[i] << " ";
//    }
//    cout << endl;
//
//    //hierarchie => model: getPhysicalGroups => physical group1,2,3,...: getEntitiesForPhysicalGroup => entity1,2,3,...: getElements => element type1,2,3,...
//
//    //looping over entities
//    for(int k = 0; k < tags.size(); ++k)
//    {
//        cout << "Entity (" << dim << "," << tags[k] << "):\n";
//        vector<int> elementTypes;
//        vector<vector<size_t>> elementTags;
//        vector<vector<size_t>> nodeTags;
//
//        model::mesh::getElements(elementTypes, elementTags, nodeTags, dim, tags[k]);
//        // looping over element types
//        for(int i = 0; i < elementTypes.size(); ++i)
//        {
//            cout << "\telement type " << elementTypes[i] << endl;
//             cout << "\t\telements: ";
//             for (int j = 0; j < elementTags[i].size(); ++j)
//             {
//                 cout << elementTags[i][j] << " ";
//             }
//             cout << endl;
//             cout << "\t\tnodes: ";
//             for (int j = 0; j < nodeTags[i].size(); ++j)
//             {
//                 cout << nodeTags[i][j] << " ";
//             }
//             cout << endl;
//        }
//    }
//
//    return 0;
//}
