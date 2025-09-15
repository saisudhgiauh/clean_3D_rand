import numpy as np
import meep as mp

def get_gaps(n_branches,offset_branches):

    '''
    inputs:
        n_branches (int) -- number of branches drawn on ridges
        offset_branches (bool) -- whether the branches are offset
    outputs:
        n_gaps (int) -- number of gaps between branches there are on the resulting structure. If branches are not offset, this is the maximum number of
            gaps on one side of the trunk only. Used for drawing the correct number of roundmasks between branches if necessary.
    '''

    if not offset_branches:
        if n_branches%2 == 0:   #even number of branches
            return (n_branches / 2) - 1
        else:       #odd number of branches
            return ((n_branches + 1) / 2) - 1
    else:
        return n_branches - 1

def branch_width_list_no_offset(taper_level,start_taper_at,max_width,n_branches):

    '''
    inputs:
        taper_level (float 0-1) -- degree of tapering
        start_taper_at (int) -- which branch in the sequence to begin the taper
        max_width (float) -- untapered branch width
        n_branches -- total number of branches on the tree
    outputs:
        widths (1D array of floats, length n_branches) -- array of the suitable branch widths, accounting for tapering
    '''
    
    widths = np.zeros(n_branches)

    if taper_level == 0:
        for i in range(n_branches):
            widths[i] = max_width
    else:
        top_width = max_width - taper_level * max_width
        for i in range(start_taper_at):     #0 to start_taper_at-1
            widths[i] = max_width
        total_gaps = get_gaps(n_branches,False)
        taper_gaps = total_gaps - get_gaps(start_taper_at,False)
        delta_w = (max_width - top_width) / taper_gaps
        for i in range(start_taper_at,n_branches):  #starts at an even number
            if i%2 == 0: #i even
                widths[i] = max_width - delta_w * ((i - (start_taper_at)) + 2) / 2
                if i+1 < n_branches:
                    widths[i+1] = max_width - delta_w * ((i - (start_taper_at)) + 2) / 2
    return widths

def branch_width_list_offset(taper_level,start_taper_at,max_width,n_branches):

    '''
    same as the function above but with offset branches
    '''

    widths = np.zeros(n_branches)

    if taper_level == 0:
        for i in range(n_branches):
            widths[i] = max_width
    else:
        top_width = max_width - taper_level * max_width
        for i in range(start_taper_at):
            widths[i] = max_width
        total_gaps = get_gaps(n_branches,True)
        taper_gaps = total_gaps - get_gaps(start_taper_at,True)
        delta_w = (max_width - top_width) / taper_gaps
        for i in range(start_taper_at,n_branches):
            widths[i] = max_width - delta_w * (i - (start_taper_at - 1))
    return widths

###a hybrid meep object created from a rectangle and a half-circle, used as a 'mask' for the
###rounding off process for branches and gaps.
###returns an array which must be +='d onto the existing geometry array, this is a bit
###cursed but meep doesn't let you create your own GeometricObject subclasses without
###editing the source so it's good enough
###the 'cent' parameter is the point at which a tangent to the half-circle is parallel to
###the far side of the block. This is usually the end of whatever you're trying to
###round off.

def horizontal_roundmask(orientation,height,width,cent,rect_mat,circ_mat):

    '''
    inputs:
        orientation (string, "left" or "right") -- the direction for the fangs of the mask to point
        height (float) -- height of mask (see diagram in the notes folder for clarification)
        width (float) -- width of mask
        cent (mp.Vector3 object) -- centre of mask
        rect_mat (mp.Medium object) -- material for the rectangular(2D)/cuboidal(3D) part of the mask
        circ_mat (mp.Medium object) -- material for the circular(2D)/cylindrical(3D) part of the mask
    outputs:   
        mask (list of GeometricObject classes) -- the MEEP objects that comprise the mask, to be appended to the simulation geometry list.
            Used for rounding off the sharp edges of branches and the spaces between branches.
    '''

    rad = height / 2
    if orientation == "right": #'fangs' point to the right
        return [
            mp.Block(
                mp.Vector3(width,height,mp.inf),
                center=(cent - mp.Vector3((width / 2) - rad,0,0)),
                material=rect_mat
            ),
            mp.Wedge(
                rad,
                height=mp.inf,
                axis=mp.Vector3(0,0,1),
                center=(cent + mp.Vector3(rad,0,0)),
                material=circ_mat,
                wedge_angle=np.pi,
                wedge_start=mp.Vector3(0,1,0)
            ) 
        ]
    elif orientation == "left": #'fangs' point to the left
        return [
            mp.Block(
                mp.Vector3(width,height,mp.inf),
                center=(cent + mp.Vector3((width / 2) - rad,0,0)),
                material=rect_mat
            ),
            mp.Wedge(
                rad,
                height=mp.inf,
                axis=mp.Vector3(0,0,1),
                center=(cent - mp.Vector3(rad,0,0)),
                material=circ_mat,
                wedge_angle=np.pi,
                wedge_start=mp.Vector3(0,-1,0)
            )
        ]
    
def vertical_roundmask(orientation,height,width,cent,rect_mat,circ_mat):

    '''
    inputs:
        orientation (string, "up" or "down") -- the direction for the fangs of the mask to point
    all other inputs and outputs same as the function above. Used for rounding off the top of the trunk.
    '''

    rad = width / 2
    if orientation == "up":     #'fangs' point upwards
        return [
            mp.Block(
                mp.Vector3(width,height,mp.inf),
                center=(cent - mp.Vector3(0,(height / 2) - rad,0)),
                material=rect_mat
            ),
            mp.Wedge(
                rad,
                height=mp.inf,
                axis=mp.Vector3(0,0,1),
                center=(cent +  mp.Vector3(0,rad,0)),
                material=circ_mat,
                wedge_angle=np.pi,
                wedge_start=mp.Vector3(-1,0,0)
            )
        ]
    elif orientation == "down":
        return [
            mp.Block(
                mp.Vector3(width,height,mp.inf),
                center=(cent + mp.Vector3(0,(height / 2) - rad,0)),
                material=rect_mat
            ),
            mp.Wedge(
                rad,
                height=mp.inf,
                axis=mp.Vector3(0,0,1),
                center=(cent - mp.Vector3(0,rad,0)),
                material=circ_mat,
                wedge_angle=np.pi,
                wedge_start=mp.Vector3(1,0,0)
            )
        ]

def roundmask_3D(height,width,length,cent,axis,fang_direction,rect_mat,cyl_mat):

    '''
    inputs:
        height (float) -- height of mask (see diagram in the notes folder for clarification)
        width (float) -- width of mask
        length (float) -- length of mask
        cent (mp.Vector3 object) -- centre of mask
        axis (mp.Vector3 object, length ignored) -- axis of the cylindrical component
        fang_direction (mp.Vector3 object, length ignored) -- direction for the fangs of the mask to point. Must be orthogonal to axis.
        rect_mat (mp.Medium object) -- material for the rectangular(2D)/cuboidal(3D) part of the mask
        cyl_mat (mp.Medium object) -- material for the circular(2D)/cylindrical(3D) part of the mask
    outputs:
        mask (list of GeometricObject classes) -- the MEEP objects that comprise the mask, to be appended to the simulation geometry list.
            Used in a similar way to horizontal_roundmask and vertical_roundmask, except this is generalised to work properly in 3D, with arbitrary
            angles and axes. Needed for drawing ridges with rotation.
    '''

    rad = height / 2
    cross = fang_direction.cross(axis)
    fang_norm = fang_direction.norm()
    return [
        mp.Block(
            mp.Vector3(height,width,length),
            e2=fang_direction,
            e3=axis,
            e1=cross,
            center=(cent + (rad - (width / 2)) * fang_direction.scale(1/fang_norm)),
            material=rect_mat,
        ),
        mp.Wedge(
            rad,
            height=length,
            axis=axis,
            center=(cent + rad * fang_direction.scale(1/fang_norm)),
            material=cyl_mat,
            wedge_angle=np.pi,
            wedge_start=-1*cross
        )
    ]