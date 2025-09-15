import numpy as np
import meep as mp
from config import *
import libs.taper_and_rounding_functions as tf

###highly painful coordinate geometry to construct the ridge dielectric
###many brain cells were lost in getting this to work as intended
###the demon of babylon disguises himself with the coat of the righteous

def get_ridge(x,eta_degrees):

    '''
    inputs:
        x (float) -- x coordinate to draw the centre of a ridge. Must be inside the cell boundaries.
        eta_degrees (float) -- angle in degrees from the vertical to draw this particular ridge at.
    outputs:
        geom_array (list of GeometricObject classes) -- the MEEP objects that comprise the ridge, to be appended to the simulation geometry list.
    '''

    geom_array = []
    
    eta = np.deg2rad(eta_degrees)
    tv = mp.Vector3(np.sin(eta),np.cos(eta),0)
    bv = mp.Vector3(np.cos(eta),-np.sin(eta),0)
    right_branch_first_centre = mp.Vector3((trunk_width / 2) + (branch_width / 2),(trunk_height / 2) - (3/2 * branch_height),0).rotate(mp.Vector3(0,0,1),-eta) + mp.Vector3(x,0,0) #this is the position of the highest branch on the trunk WITHOUT TAPERING, subsequent branch positions are generated from this
    left_branch_first_centre_offset = right_branch_first_centre - mp.Vector3((branch_width + trunk_width),(branch_height / 2) + (branch_sep / 2),0).rotate(mp.Vector3(0,0,1),-eta)
    left_branch_first_centre_no_offset = right_branch_first_centre - mp.Vector3(branch_width + trunk_width,0,0).rotate(mp.Vector3(0,0,1),-eta)

        #round centers
        #note: the seemingly random "+ 0.01 * branch_sep" bits on the height of each mask are to
        #prevents strange artefacts due to imperfect tilings of geometric components.
        #can be removed without changing the spectrum too much if this starts causing problems

    if round_centres == True:
        right_gap_first_mask = right_branch_first_centre - mp.Vector3(branch_width / 2,(branch_height / 2) + (branch_sep / 2),0).rotate(mp.Vector3(0,0,1),-eta)
        left_gap_first_mask_offset = left_branch_first_centre_offset - mp.Vector3(-1 * branch_width / 2,(branch_height / 2) + (branch_sep / 2),0).rotate(mp.Vector3(0,0,1),-eta)
        left_gap_first_mask_no_offset = left_branch_first_centre_no_offset - mp.Vector3(-1 * branch_width / 2,(branch_height / 2) + (branch_sep / 2),0).rotate(mp.Vector3(0,0,1),-eta)
        if offset_branches:
            n_gaps = int(tf.get_gaps(n_branches,offset_branches))
            first_left = left_gap_first_mask_offset
        else:
            n_gaps = 2 * int(tf.get_gaps(n_branches,offset_branches))
            first_left = left_gap_first_mask_no_offset
        for i in range(1,n_gaps+1):
            if i%2 != 0:   #i odd, right gaps
                ith_mask_cent = right_gap_first_mask - (i-1)/2 * mp.Vector3(0,branch_sep + branch_height,0).rotate(mp.Vector3(0,0,1),-eta)
                geom_array += tf.roundmask_3D(
                    height=branch_sep * 1.01,
                    width=branch_sep / 2,
                    length=mp.inf,
                    cent=ith_mask_cent,
                    axis=mp.Vector3(0,0,1),
                    fang_direction=bv,
                    rect_mat=wing_mat,
                    cyl_mat=air
                )
            else:   #i even, left gaps
                ith_mask_cent = first_left - (i-2)/2 * mp.Vector3(0,branch_sep + branch_height,0).rotate(mp.Vector3(0,0,1),-eta)
                geom_array += tf.roundmask_3D(
                    height=branch_sep * 1.01,
                    width=branch_sep / 2,
                    length=mp.inf,
                    cent=ith_mask_cent,
                    axis=mp.Vector3(0,0,1),
                    fang_direction=-1*bv,
                    rect_mat=wing_mat,
                    cyl_mat=air
                )    

        #now we add the branches by appending additional meep objects to the geom_array list
        #note: at the moment this is hard to read. should change this in future - there are
        #a few common statements that could be aliased, and some bits are repeated unnecessarily.

    if n_branches != 0:
        if offset_branches:
            branch_widths_list = tf.branch_width_list_offset(taper_level,start_taper_at,branch_width,n_branches)
            first_left = left_branch_first_centre_offset
        else:
            branch_widths_list = tf.branch_width_list_no_offset(taper_level,start_taper_at,branch_width,n_branches)
            first_left = left_branch_first_centre_no_offset
        for i in range(1,n_branches+1):
            if i%2 != 0:     #i odd, right branches
                ith_branch_cent = right_branch_first_centre - (i-1)/2 * mp.Vector3(0,(branch_sep + branch_height),0).rotate(mp.Vector3(0,0,1),-eta) - mp.Vector3((branch_width - np.flip(branch_widths_list)[i-1]) / 2,0,0).rotate(mp.Vector3(0,0,1),-eta)
                geom_array.append(
                    mp.Block(
                        mp.Vector3(np.flip(branch_widths_list)[i-1],branch_height,mp.inf),
                        center=ith_branch_cent,
                        e1=bv,
                        e2=tv,
                        material=wing_mat
                    )
                )
                if round_edges:
                    geom_array += tf.roundmask_3D(
                        height=branch_height,
                        width=branch_height/2 * 1.01,
                        length=mp.inf,
                        cent=ith_branch_cent + mp.Vector3(np.flip(branch_widths_list)[i-1] / 2,0,0).rotate(mp.Vector3(0,0,1),-eta),
                        axis=mp.Vector3(0,0,1),
                        fang_direction=-1*bv,
                        rect_mat=air,
                        cyl_mat=wing_mat
                    )
            else:   #i even, left branches
                ith_branch_cent = first_left - (i-2)/2 * mp.Vector3(0,(branch_sep + branch_height),0).rotate(mp.Vector3(0,0,1),-eta) + mp.Vector3((branch_width - np.flip(branch_widths_list)[i-1]) / 2,0,0).rotate(mp.Vector3(0,0,1),-eta)
                geom_array.append(
                    mp.Block(
                        mp.Vector3(np.flip(branch_widths_list)[i-1],branch_height,mp.inf),
                        center=ith_branch_cent,
                        e1=bv,
                        e2=tv,
                        material=wing_mat
                    )
                )
                if round_edges:
                    geom_array += tf.roundmask_3D(
                        height=branch_height,
                        width=branch_height/2 * 1.01,
                        length=mp.inf,
                        cent=ith_branch_cent - mp.Vector3(np.flip(branch_widths_list)[i-1] / 2,0,0).rotate(mp.Vector3(0,0,1),-eta),
                        axis=mp.Vector3(0,0,1),
                        fang_direction=bv,
                        rect_mat=air,
                        cyl_mat=wing_mat
                    )

        #add basal layer, trunk and overlap

    d1 = (cell_z / 2) * np.tan(np.deg2rad(overlap_angle))
    d2 = (cell_z / 2 + overlap) * np.tan(np.deg2rad(overlap_angle))
    vertices1 = [   #vertices of the triangular overlap prism
        mp.Vector3(-trunk_width / 2,trunk_height / 2,0).rotate(mp.Vector3(0,0,1),-eta) + mp.Vector3(x,0,0),
        mp.Vector3(-trunk_width / 2,trunk_height / 2,-cell_z / 2).rotate(mp.Vector3(0,0,1),-eta) + mp.Vector3(x,0,0),
        mp.Vector3(-trunk_width / 2,trunk_height / 2 + d1,-cell_z / 2).rotate(mp.Vector3(0,0,1),-eta) + mp.Vector3(x,0,0)
    ]
    vertices2 = [   #vertices of the quad overlap prism
        mp.Vector3(-trunk_width / 2,trunk_height / 2,0).rotate(mp.Vector3(0,0,1),-eta) + mp.Vector3(x,0,0),
        mp.Vector3(-trunk_width / 2,trunk_height / 2,cell_z / 2).rotate(mp.Vector3(0,0,1),-eta) + mp.Vector3(x,0,0),
        mp.Vector3(-trunk_width / 2,trunk_height / 2 + d1,cell_z / 2).rotate(mp.Vector3(0,0,1),-eta) + mp.Vector3(x,0,0),
        mp.Vector3(-trunk_width / 2,trunk_height / 2 + d1 + d2,-overlap).rotate(mp.Vector3(0,0,1),-eta) + mp.Vector3(x,0,0)
    ]

    geom_array += [
        mp.Block(       #base
            mp.Vector3(mp.inf,basal_height,mp.inf),
            center=mp.Vector3(0,-1*((cell_y / 2) - (basal_height / 2)),0),
            material=wing_mat
        ),
        mp.Block(       #main trunk
            mp.Vector3(trunk_width,trunk_height,mp.inf),
            center=mp.Vector3() + mp.Vector3(x,0,0),
            e1=bv,
            e2=tv,
            material=wing_mat
        )]

    geom_array += [mp.Prism(       #overlap triangular prism
            vertices1,
            height=trunk_width,
            axis=bv,
            material=wing_mat
        )]

    if round_top:      #round top of triangular section
        geom_array += tf.roundmask_3D(
            height=trunk_width,
            width=(3/2) * trunk_width,
            length=(vertices1[0] - vertices1[2]).norm() * 1.01,
            cent=(1/2) * (vertices1[2] + vertices1[0]) + mp.Vector3(trunk_width / 2,0,0).rotate(mp.Vector3(0,0,1),-eta),
            axis=mp.Vector3(0,-np.sin(np.deg2rad(overlap_angle)),np.cos(np.deg2rad(overlap_angle))).rotate(mp.Vector3(0,0,1),-eta),
            fang_direction=mp.Vector3(0,-np.cos(np.deg2rad(overlap_angle)),-np.sin(np.deg2rad(overlap_angle))).rotate(mp.Vector3(0,0,1),-eta),
            rect_mat=air,
            cyl_mat=wing_mat
        )

    geom_array += [mp.Prism(       #overlap quad prism
            vertices2,
            height=trunk_width,
            axis=bv,
            material=wing_mat
        )]

    if round_top:    
        if overlap != 0:
            overhang_angle = np.arctan((d1+d2)/overlap)
        else:
            overhang_angle = 0
        geom_array += tf.roundmask_3D(     #round overhang of quad prism
            height=trunk_width,
            width=(1/2) * trunk_width,
            length=(vertices1[0] - vertices2[3]).norm(),
            cent=(1/2) * (vertices1[0] + vertices2[3]) + mp.Vector3(trunk_width / 2,0,0).rotate(mp.Vector3(0,0,1),-eta),
            axis=mp.Vector3(0,np.sin(overhang_angle),-np.cos(overhang_angle)).rotate(mp.Vector3(0,0,1),-eta),
            fang_direction=mp.Vector3(0,np.cos(overhang_angle),np.sin(overhang_angle)).rotate(mp.Vector3(0,0,1),-eta),
            rect_mat=air,
            cyl_mat=wing_mat
        )
        geom_array += tf.roundmask_3D(     #round top of quad prism
            height=trunk_width,
            width=(3/2) * trunk_width,      
            length=(vertices2[2] - vertices2[3]).norm(),
            cent=(1/2) * (vertices2[2] + vertices2[3]) + mp.Vector3(trunk_width / 2,0,0).rotate(mp.Vector3(0,0,1),-eta),
            axis=mp.Vector3(0,-np.sin(np.deg2rad(overlap_angle)),np.cos(np.deg2rad(overlap_angle))).rotate(mp.Vector3(0,0,1),-eta),
            fang_direction=mp.Vector3(0,-np.cos(np.deg2rad(overlap_angle)),-np.sin(np.deg2rad(overlap_angle))).rotate(mp.Vector3(0,0,1),-eta),
            rect_mat=air,
            cyl_mat=wing_mat
        )

    geom_array += [mp.Block(       #fix bottom of trunk
            mp.Vector3(trunk_width,basal_height/2,mp.inf),
            center=mp.Vector3((-trunk_height/2)*np.sin(eta),(-trunk_height/2)*np.cos(eta),0) + mp.Vector3(x,0,0),
            e1=bv,
            e2=tv,
            material=wing_mat
        )]
    
    return geom_array