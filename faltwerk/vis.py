from itertools import product
# import pdb
# import warnings

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
# plt.get_cmap('viridis')

import py3Dmol


def plot_alphafold(af, width=400, height=300):

    styles = {
        'A': ('black', 1),
        'B': ('grey', 0.50),
        'C': ('grey', 0.50),
        'D': ('grey', 0.50),
        'E': ('grey', 0.50),
    }
    view=py3Dmol.view(width=width, height=height)
    for m in af.models.values():
        view.addModel(m.to_stream(),'pdb')

    view.setBackgroundColor('white')

    for k, (color, opacity) in styles.items():
        view.setStyle({'chain': k},{'cartoon': {'color': color, 'opacity': opacity}})

    view.zoomTo()
    return view


def map_colors(v, palette='viridis', limits=None, offset=0):
    '''
    offset .. in 3Dmol.js indexing of eg residues starts at 1, so we need to
    start indexing at 1 in the corresponding color map, ie offset=1
    '''
    if not limits:
        norm = matplotlib.colors.Normalize(vmin=min(v), vmax=max(v), clip=True)
    else:
        assert len(limits) == 2, 'Please supply limits like [0, 1], exit.'
        norm = matplotlib.colors.Normalize(
            vmin=limits[0], vmax=limits[1], clip=True)
    
    mapper = cm.ScalarMappable(norm=norm, cmap=plt.get_cmap(palette))
    
    colors = [matplotlib.colors.to_hex(mapper.to_rgba(i)) for i in v]
    cmap = {k+offset:v for k, v in enumerate(colors)}
    return cmap


def plot_annotation(model, label, palette='viridis', surface=False, opacity=1., width=400, height=300, limits=None):
    '''
    Available color maps:

    https://matplotlib.org/stable/tutorials/colors/colormaps.html
    '''
    assert label in model.annotation.keys(), 'Label not found in annotation'
    stream = model.to_stream()
    
    # https://stackoverflow.com/questions/28752727/map-values-to-colors-in-matplotlib
    anno = model.annotation[label]
    
    d = map_colors(anno, palette, limits)

    # mn = min(anno)
    # mx = max(anno)
    
    # norm = matplotlib.colors.Normalize(vmin=mn, vmax=mx, clip=True)
    # mapper = cm.ScalarMappable(norm=norm, cmap=plt.get_cmap(palette))
    
    # # Map colors to residues
    # d = {k+1: matplotlib.colors.to_hex(mapper.to_rgba(v)) for k, v in enumerate(anno)}
    
    view = py3Dmol.view(width=width, height=height)
    # view.addModelsAsFrames(stream)
    view.addModel(stream, 'pdb')

    i = 0
    l = []
    for line in stream.split("\n"):
        split = line.split()
        if len(split) == 0 or split[0] != "ATOM":
            continue

        # print(split)
        color = d[int(split[5])]
        l.append(color)
        # idx = int(split[1])
        view.setStyle({'model': -1, 'serial': i+1}, {'cartoon': {'color': color}})
        i += 1
        
    if surface:
        # map_ = {(i+1): j for i, j in zip(range(len(model)), map_colors(anno, palette, limits))}
        view.addSurface(py3Dmol.VDW, {'opacity': opacity, 'colorscheme': {'prop': 'resi', 'map': d}})

    view.zoomTo()
    return view


def plot_superposition(models, colors, width=400, height=300):
    
    view=py3Dmol.view(width=width, height=height)
    for m in models:
        view.addModel(m.to_stream(), 'pdb')

    view.setBackgroundColor('white')

    for chain, color in colors.items():
        view.setStyle({'chain': chain}, {'cartoon': {'color': color}})

    view.zoomTo()
    return view


class Layout():
    '''
    We should be able to overwrite say ribbon with then sphere and lastly surface
    select - style - select - style ... layer on top
    
    https://stackoverflow.com/questions/41817578/basic-method-chaining
    '''
    def __init__(self, fold=None, panel_size=(400, 300), grid=(1, 1), linked=True, src='https://3dmol.org/build/3Dmol.js'):
        self.fold = fold
        
        # Create panel map, which tracks which vis is applied in which panel of the layout
        rows = [i for i in range(grid[0])]
        cols = [i for i in range(grid[1])]
        self.build = {p: [] for p in product(rows, cols)}
        
        # "panel", not "facet" -- https://ggplot2.tidyverse.org/reference/facet_grid.html
        view = py3Dmol.view(
            width=panel_size[0]*grid[1],
            height=panel_size[1]*grid[0],
            linked=linked,
            viewergrid=grid,
            js=src)
        stream = self.fold.to_stream()
        view.addModel(stream, 'pdb')
        self.view = view

    def select(self, residues=[], elements=[], chain=''):
        '''
        - fold.structure > model > chain > residue > atom
        There are three items one might want to regularly select: chains, residues and atoms.
        
        residues .. [1, 2, 89, 90, 91]
        elements .. ['CA', 'N', 'O']
        '''
        if residues:
            selection = {'resi': [str(i) for i in residues]}
        else:
            selection = {'resi': [f'1-{len(self.fold)}']}
            
        if elements:
            selection.update({'atom': elements})
        
        if chain:
            assert chain in [i.id for i in self.fold.structure.get_chains()]
            selection.update({'chain': chain})
        
        return selection
        
    def geom_ribbon(self, key=None, selection=None, color='grey', palette='viridis', color_map_limits=None, thickness=0.4, opacity=1, ribbon=False, arrows=False, tubes=False, style='rectangle', panel=(0, 0)):
        '''
        - https://3dmol.csb.pitt.edu/doc/types.html#ViewerGridSpec

        https://github.com/matplotlib/matplotlib/blob/f6e0ee49c598f59c6e6cf4eefe473e4dc634a58a/examples/color/colormap_reference.py#L8
        
        A reversed version of each of these colormaps is available by appending
        _r to the name, e.g., viridis_r.

        '''        
        if not selection:
            selection = self.select()

        # Properties directly passed on to 3Dmol.js:
        # - https://3dmol.csb.pitt.edu/doc/types.html#CartoonStyleSpec
        # - https://3dmol.csb.pitt.edu/doc/types.html#CrossStyleSpec
        assert style in ['trace', 'oval', 'rectangle', 'parabola', 'edged']
        properties = {
            'thickness': thickness,
            'opacity': opacity,
            'ribbon': ribbon,
            'arrows': arrows,
            'tubes': tubes,
            'style': style,
        }
        
        if not key:
            style = {'cartoon': {'color': color, **properties}}
            # self.build[panel].append([selection, style])
            _ = self.add_to_build('ribbon', panel, [selection, style])
            return self
      
        anno = self.fold.annotation[key]
        # It could be just one value, in which case no palette and color mapping required
        n_colors = len(set(anno))
        
        if n_colors > 1:
            cmap = map_colors(anno, palette, color_map_limits, offset=1)
            # Maps positions to colors
            # ... 668: '#ff2211', 669: '#ff1f10', 670: '#ff1f10' ...
        
        # 'color': '#ff2211', 'color': 'green'
        style = {'cartoon': {'opacity': opacity, 'colorscheme': {'prop': 'resi', 'map': cmap}}}
        # _ = self.add_to_build(panel, selection, style)
        _ = self.add_to_build('ribbon', panel, [selection, style])
        # print(selection)
        
        #for i in self.fold.structure.get_atoms():
        #    # We map colors to residues but then color every atom (C, N, ...)
        #    ix = i.get_parent().id[1]    
        #    try:
        #        color =  {'color': cmap[ix]}
        #    except NameError:
        #        color = 'grey'
        #        
        #    style = {'cartoon': color | properties}
        #    selection = {'model': -1, 'serial': i.get_serial_number()}
        #    _ = self.add_to_build(panel, selection, style)   
        return self
    
    # TODO: this is only half done (pass to build dict, render() later etc.)
    # TODO: can we select only parts for surface vis? think so (see chain C selection):
    # https://3dmol.csb.pitt.edu/doc/types.html#SurfaceStyleSpec
    # viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.85,voldata: data, color:'red'},{chain:'C'});
    # def geom_ribbon(self, key=None, selection=None, color='grey', palette='viridis', color_map_limits=None, thickness=0.4, opacity=1, ribbon=False, arrows=False, tubes=False, style='rectangle', panel=(0, 0)):
    def geom_surface(self, key=None, selection=None, surface_type='MS', opacity=1, color='white', palette='viridis', color_map_limits=None, panel=(0, 0)):
        '''
        https://3dmol.csb.pitt.edu/doc/types.html#SurfaceStyleSpec
        
        SAS .. solvent access?
        VWD .. van der Waals ?
        SES .. electrostatic?
        MS  .. no idea
        '''
        assert surface_type in ['SAS', 'VWD', 'SES', 'MS']
        
        properties = {
            'opacity': opacity,
        }
        
        if not selection:
            selection = self.select()

        if not key:
            # self.view.addSurface('VDW', {'opacity': opacity, 'color': 'blue'}, selection, viewer=panel)
            style = {'color': color, **properties}
            _ = self.add_to_build('surface', panel, [surface_type, style, selection])
            # TODO
            # _ = self.add_to_build('surface', panel, ['VWD', style, selection])
            # self.build
            return self
      
        anno = self.fold.annotation[key]  
        # anno = self.fold.annotation['position']
        #cmap = map_colors(anno, 'viridis', None)
        
        # cmap = {k-1: v for k, v in cmap.items()}
        # print(cmap)
        cmap = map_colors(anno, palette, color_map_limits, offset=1)
        # 'colorscheme': 'yellowCarbon'
        # 'colorscheme': {'prop': 'resi', 'map': cmap}
        # self.view.addSurface(py3Dmol.VDW, {'opacity': opacity, 'colorscheme': {'prop': 'resi', 'map': cmap}}, viewer=panel)
        
        # Note how the selection here goes _after_ the style spec; not sure if this is by design or accidental.
        # self.view.addSurface('VDW', {'opacity': opacity, 'colorscheme': {'prop': 'resi', 'map': cmap}}, selection, viewer=panel)
        style = {'colorscheme': {'prop': 'resi', 'map': cmap}, **properties}
        _ = self.add_to_build('surface', panel, [surface_type, style, selection])
        # pdb.set_trace()
        return self
    
    def geom_sphere(self, selection=None, key='position', color='black', radius=2, panel=(0, 0)):
        
        if not selection:
            selection = self.select()
    
        properties = {
            'radius': radius
        }
    
        style = {'sphere': {'color': color, **properties}}
        _ = self.add_to_build('sphere', panel, [selection, style])
        return self
    
    def add_to_build(self, method, panel, spec):
        '''
        We go through this method instead of having eg geom_ribbon add its spec directly to the vis
        build in order to have a single checkpoint.
        '''
        assert method in ['ribbon', 'surface', 'sphere']
        try:
            # if panel provided as list, cast to tuple,
            # otherwise TypeError: unhashable type: 'list'
            self.build[tuple(panel)].append([method, spec])
        except KeyError:
            raise KeyError(f'{panel} is not a valid panel position')
        return None
    
    # TODO
    def label():
        self._label = True
        pass
    
    # def __call__(self):
    #     print('Wha?')
    #     return None
    
    def __getitem__(self, key):
        return self.build[key]
    
    def render(self):
        '''
        Build visualisation.
        '''
        # print(self.build)
        for panel, v in self.build.items():
            for method, spec in v:
                if method in ['ribbon', 'sphere']:
                    self.view.setStyle(*spec, viewer=panel)
                elif method == 'surface':
                    self.view.addSurface(*spec, viewer=panel)
                else:
                    raise ValueError(f'Style type "{method}" not implemented')
        
        self.view.zoomTo()
        return self.view
