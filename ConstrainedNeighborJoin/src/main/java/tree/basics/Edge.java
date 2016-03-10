/* ***** BEGIN LICENSE BLOCK *****
 *
 * Copyright (c) 2005-2007 Universidade de Sao Paulo, Sao Carlos/SP, Brazil.
 * All Rights Reserved.
 *
 * This file is part of Projection Explorer (PEx).
 *
 * How to cite this work:
 *  
@inproceedings{paulovich2007pex,
author = {Fernando V. Paulovich and Maria Cristina F. Oliveira and Rosane 
Minghim},
title = {The Projection Explorer: A Flexible Tool for Projection-based 
Multidimensional Visualization},
booktitle = {SIBGRAPI '07: Proceedings of the XX Brazilian Symposium on 
Computer Graphics and Image Processing (SIBGRAPI 2007)},
year = {2007},
isbn = {0-7695-2996-8},
pages = {27--34},
doi = {http://dx.doi.org/10.1109/SIBGRAPI.2007.39},
publisher = {IEEE Computer Society},
address = {Washington, DC, USA},
}
 *  
 * PEx is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free 
 * Software Foundation, either version 3 of the License, or (at your option) 
 * any later version.
 *
 * PEx is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 *
 * This code was developed by members of Computer Graphics and Image
 * Processing Group (http://www.lcad.icmc.usp.br) at Instituto de Ciencias
 * Matematicas e de Computacao - ICMC - (http://www.icmc.usp.br) of 
 * Universidade de Sao Paulo, Sao Carlos/SP, Brazil. The initial developer 
 * of the original code is Fernando Vieira Paulovich <fpaulovich@gmail.com>.
 *
 * Contributor(s): Rosane Minghim <rminghim@icmc.usp.br>
 *
 * You should have received a copy of the GNU General Public License along 
 * with PEx. If not, see <http://www.gnu.org/licenses/>.
 *
 * ***** END LICENSE BLOCK ***** */

package tree.basics;

import java.awt.Color;
import java.io.Serializable;

/**
 * @author Fernando Vieira Paulovich
 *
 * This class represents a edge on the map.
 */
public class Edge implements Comparable, Serializable {

    public static final float NO_SIZE = -1;
    private static final long serialVersionUID = 1L;

    /**
     * Constructor of the edge
     *      
     * @param source The first vertex
     * @param target The second vertex
     * @param weight The edge's weight
     */
    public Edge(int source, int target, float weight) {
        this(source, target);
        this.weight = weight;
    }

    /**
     * Constructor of the edge
     * 
     * @param source The first vertex
     * @param target The second vertex
     */
    public Edge(int source, int target) {
        this.source = source;
        this.target = target;
    }

    /**
     * Return the color of the edge
     * @return The color of the edge
     */
    public Color getColor() {
        return this.color;
    }

    /**
     * Changes the color of the edge
     * @param color
     */
    public void setColor(Color color) {
        this.color = color;
    }

    /**
     * Return the source of the edge
     * @return The source of the edge
     */
    public int getSource() {
        return this.source;
    }

    /**
     * Return the target of the edge
     * @return The target of the edge
     */
    public int getTarget() {
        return this.target;
    }

    /**
     * Return the weight of the edge
     * @return The weight of the edge
     */
    public float getWeight() {
        return weight;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof Edge) {
            Edge e = (Edge) obj;
            return (((this.source == e.source) && (this.target == e.target)) ||
                    ((this.source == e.target) && (this.target == e.source)));
        }
        return false;
    }

    @Override
    public int hashCode() {
        int hash = 3 + 5 * this.source;
        hash += 7 * this.target;
        return hash;
    }

    @Override
    public int compareTo(Object o) {
        long source_aux = 0;
        long target_aux = 0;

        if (this.source < this.target) {
            source_aux = this.source;
            target_aux = this.target;
        } else {
            source_aux = this.target;
            target_aux = this.source;
        }

        long sourceComp = 0;
        long targetComp = 0;
        if (((Edge) o).source < ((Edge) o).target) {
            sourceComp = ((Edge) o).source;
            targetComp = ((Edge) o).target;
        } else {
            sourceComp = ((Edge) o).target;
            targetComp = ((Edge) o).source;
        }

        if (source_aux - sourceComp < 0) {
            return -1;
        } else if (source_aux - sourceComp > 0) {
            return 1;
        } else {
            if (target_aux - targetComp < 0) {
                return -1;
            } else if (target_aux - targetComp > 0) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    private float weight = Edge.NO_SIZE;
    private Color color = Color.WHITE; //Color of the edge
    private int source; //The first vertex of the edge
    private int target; //The second vertex of the edge
}
