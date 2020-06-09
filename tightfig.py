# Generated with SMOP  0.41-beta
from libsmop import *
# ddd/tightfig.m

    
@function
def tightfig(hfig=None,*args,**kwargs):
    varargin = tightfig.varargin
    nargin = tightfig.nargin

    # tightfig: Alters a figure so that it has the minimum size necessary to
# enclose all axes in the figure without excess space around them.
# 
# Note that tightfig will expand the figure to completely encompass all
# axes if necessary. If any 3D axes are present which have been zoomed,
# tightfig will produce an error, as these cannot easily be dealt with.
# 
# hfig - handle to figure, if not supplied, the current figure will be used
# instead.
    
    if nargin == 0:
        hfig=copy(gcf)
# ddd/tightfig.m:13
    
    # There can be an issue with tightfig when the user has been modifying
    # the contnts manually, the code below is an attempt to resolve this,
    # but it has not yet been satisfactorily fixed
#     origwindowstyle = get(hfig, 'WindowStyle');
    set(hfig,'WindowStyle','normal')
    
    # get all the axes handles note this will also fetch legends and
    # colorbars as well
    hax=findall(hfig,'type','axes')
# ddd/tightfig.m:26
    
    # later
    origaxunits=get(hax,'Units')
# ddd/tightfig.m:30
    
    set(hax,'Units','centimeters')
    
    if numel(hax) > 1:
        #         fsize = cell2mat(get(hax, 'FontSize'));
        ti=cell2mat(get(hax,'TightInset'))
# ddd/tightfig.m:38
        pos=cell2mat(get(hax,'Position'))
# ddd/tightfig.m:39
    else:
        #         fsize = get(hax, 'FontSize');
        ti=get(hax,'TightInset')
# ddd/tightfig.m:42
        pos=get(hax,'Position')
# ddd/tightfig.m:43
    
    
    # ensure very tiny border so outer box always appears
    ti[ti < 1.25]=1.25
# ddd/tightfig.m:47
    
    # they are not being viewed in any of the 2d directions
    views2d=concat([[0,90],[0,0],[90,0]])
# ddd/tightfig.m:51
    for i in arange(1,numel(hax)).reshape(-1):
        set(hax(i),'LooseInset',ti(i,arange()))
        #         set(hax(i), 'LooseInset', [0,0,0,0]);
        # get the current viewing angle of the axes
        az,el=view(hax(i),nargout=2)
# ddd/tightfig.m:59
        iszoomed=strcmp(get(hax(i),'CameraViewAngleMode'),'manual')
# ddd/tightfig.m:62
        is2d=all(bsxfun(eq,concat([az,el]),views2d),2)
# ddd/tightfig.m:65
        if iszoomed and logical_not(any(is2d)):
            error('TIGHTFIG:haszoomed3d','Cannot make figures containing zoomed 3D axes tight.')
    
    
    # we will move all the axes down and to the left by the amount
    # necessary to just show the bottom and leftmost axes and labels etc.
    moveleft=min(pos(arange(),1) - ti(arange(),1))
# ddd/tightfig.m:75
    movedown=min(pos(arange(),2) - ti(arange(),2))
# ddd/tightfig.m:77
    
    # encompass the topmost and rightmost axes and lables
    figwidth=max(pos(arange(),1) + pos(arange(),3) + ti(arange(),3) - moveleft)
# ddd/tightfig.m:81
    figheight=max(pos(arange(),2) + pos(arange(),4) + ti(arange(),4) - movedown)
# ddd/tightfig.m:83
    
    for i in arange(1,numel(hax)).reshape(-1):
        set(hax(i),'Position',concat([pos(i,arange(1,2)) - concat([moveleft,movedown]),pos(i,arange(3,4))]))
    
    
    origfigunits=get(hfig,'Units')
# ddd/tightfig.m:92
    set(hfig,'Units','centimeters')
    
    figpos=get(hfig,'Position')
# ddd/tightfig.m:97
    set(hfig,'Position',concat([figpos(1),figpos(2),figwidth,figheight]))
    
    set(hfig,'PaperUnits','centimeters')
    set(hfig,'PaperSize',concat([figwidth,figheight]))
    set(hfig,'PaperPositionMode','manual')
    set(hfig,'PaperPosition',concat([0,0,figwidth,figheight]))
    
    if logical_not(iscell(origaxunits)):
        origaxunits=cellarray([origaxunits])
# ddd/tightfig.m:109
    
    for i in arange(1,numel(hax)).reshape(-1):
        set(hax(i),'Units',origaxunits[i])
    
    set(hfig,'Units',origfigunits)
    #      set(hfig, 'WindowStyle', origwindowstyle);
    
    return hfig
    
if __name__ == '__main__':
    pass
    