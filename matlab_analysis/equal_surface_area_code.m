%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Voronoi tessellation with equilibration Code
%
% Step 1 - Generate initial centroids (either randomly or
%          using a predefined method, N centroids = N parcels)
%
% Step 2 - Generate Voronoi tessellation using those centroids
%
% Step 3 - Jiggle centroids around until they are equidistant from each
%          other. This will typically gives an equal surface area
%          parcellation.
%          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variable Directory
% number_centroid - number of centroids
% data_dir - directory for mesh
% verts - vertices of mesh
% faces - faces of mesh
% mask - mask (0 for outside mask, 1 for inside mask)
% edges - edges of mask used in main code
% dist - matrix containing the distances between vertices on the mesh
% NNs - list containing the nearest neighbour information about the mesh
% initial_centroids - inital guess for the location of the centroids
%                     (should be fixed for consistent results)
% new_parc - generated parcellation on full mesh (inc. masked region)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Example use of the code -  Generating 232 parcels for a human adult lh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define number of centroids
number_centroids=100; 

% Import mesh and mask
% data_dir=[pwd,'\fsLR_32k_midthickness-lh.vtk'];
% [verts,faces]=read_vtk(data_dir,[]);

G = gifti('./dHCP40wk_surface/week-40_hemi-right_space-dhcpSym_dens-32k_wm.surf.gii');

verts = double(G.vertices);
faces = double(G.faces);

MedMask = gifti('./dHCP40wk_surface/week-40_hemi-left_space-dhcpSym_dens-32k_desc-medialwallsymm_mask.shape.gii');

mask=double(MedMask.cdata);

% Dalculate mesh metrics used in main code (this section taks the longest
% to run)
if exist('mask') == 1
    edges=edges_wb(verts,faces,mask);
    [verts,faces]=wall_cut(verts,faces,find(mask==1));
else
    edges=[];
end
dist=dist_graph(verts,faces);
dist=1/2*(dist+dist');
NNs=neighs(verts,faces);

% Generate initial seed for centroids 
initial_centroids=[];
removed=[];
for i=1:number_centroids
    brek=0;
    while brek==0
        ii=randi(max(size(verts)),1);
        if ~ismember(ii,removed)
            initial_centroids(end+1)=ii;
            removed(end+1:end+numel(NNs{ii})+1)=[NNs{ii};ii];
            brek=1;
        end
    end
end
clear removed brek

% Run main code
new_centroids=mesh_subsegment(verts,initial_centroids,edges,NNs,dist);
[new_parc]=startpts_to_parcer(new_centroids,NNs,dist,edges,mask);

function [full_parc]=startpts_to_parcer(start_points,NNs,dist,edges,mask)
    [up_parc,~]=perform_badest_updater(start_points,NNs,dist,edges);
    uQ=start_points;
    parc=zeros(size(up_parc));
    for i=1:numel(uQ)
        parc(up_parc==uQ(i))=i;
    end
    full_parc=zeros(max(size(mask)),1);
    full_parc(mask==1)=parc;
end

function [start_points]=mesh_subsegment(verts,start_points,edges,NNs,dist)
    if size(verts,1)<size(verts,2)
        verts = verts';
    end

    dsold=0;
    counter=0;
    dslist=[];
    for kk=1:5000
        start_points2=start_points;
        [parc,D]=perform_badest_updater(start_points,NNs,dist,edges);
        %start_points=[];
        uQ=start_points;
        for j=1:numel(uQ)
            pts=find(parc==uQ(j));
            [~,k] = max(D(pts));
            start_points(j) = pts(k);
        end
        ds=sum(sqrt(sum((verts(start_points2(:),:)-verts(start_points(:),:)).^2,2)));
        if abs(ds-dsold)<1
            counter=counter+1;
        else
            counter=0;
        end
        if counter==15
            break
        end
        dsold=ds;
        if kk>15
            dslist(end+1)=ds;
            [~,bb,~]=mode(dslist);
            if bb>15
                break
            end
        end
    end
end
        

function [verts,facs]=wall_cut(verts,facs,mask_list)
    mask_list=mask_list(:);
    if size(verts,1)<size(verts,2)
        verts = verts';
    end
    if size(facs,1)<size(facs,2)
        facs = facs';
    end
    NNs=neighs(verts,facs);
    new_mask_list=[];
    for i=1:size(mask_list,1)
        aaa=NNs{mask_list(i)};
        vals=double(ismember(aaa,mask_list));
        if sum(vals)>1
            new_mask_list(end+1)=mask_list(i);
        end
    end

    indy=zeros(size(verts,1),1);
    indy(new_mask_list)=new_mask_list;
    aa=find(indy~=0);
    
    verts2=verts(new_mask_list,:);
    facs2=[];
    for i=1:size(facs,1)
        if sum(double(ismember(facs(i,:),new_mask_list)))==3
            gg=facs(i,:);
            fff=mod(find(aa==double(facs(i,:)))-1,size(aa,1))+1;
            facs2(end+1,:)=fff;
        end
    end
    facs=facs2;
    verts=verts2;
end

function [dist,G]=dist_graph(verts,facs)
    if size(verts,1)<size(verts,2)
        verts = verts';
    end
    if size(facs,1)<size(facs,2)
        facs = facs';
    end

    stt=size(facs,1);
    ss=[];
    tt=[];
    weights=[];
    for i=1:size(verts,1)
        aaa=unique(facs(mod(find(facs==i)-1,stt)+1,:));
        aaa=aaa(aaa~=i);
        for j=1:numel(aaa)
            ss(end+1)=i;
            tt(end+1)=aaa(j);
            weights(end+1)=sqrt((verts(i,1)-verts(aaa(j),1)).^2+(verts(i,2)-verts(aaa(j),2)).^2+(verts(i,3)-verts(aaa(j),3)).^2);
        end
    end
    
    G=graph(ss,tt,weights);
    dist = distances(G);

end

function NNs=neighs(verts,facs)
    if size(verts,1)<size(verts,2)
        verts = verts';
    end
    if size(facs,1)<size(facs,2)
        facs = facs';
    end
    stt=size(facs,1);
    NNs={};
    for i=1:size(verts,1)
        aaa=unique(facs(mod(find(facs==i)-1,stt)+1,:));
        aaa=aaa(aaa~=i);
        NNs{i}=aaa;
    end
end

function edges=edges_wb(verts_wb,facs_wb,mask_list)
    if size(verts_wb,1)<size(verts_wb,2)
        verts_wb = verts_wb';
    end
    if size(facs_wb,1)<size(facs_wb,2)
        facs_wb = facs_wb';
    end

    stt=size(facs_wb,1);
    edges=zeros(size(mask_list));
    for i=1:size(verts_wb,1)
        aaa=unique(facs_wb(mod(find(facs_wb==i)-1,stt)+1,:));
        ddd=unique(mask_list(aaa));
        if numel(ddd)>1
        edges(i)=1;
        end
    end
    edges=edges(mask_list==1);
    edges=find(edges);
end

 function [assign,D]=perform_badest_updater(start_points,NNs,dist,edges)
    dists=dist(start_points,:);
    map=zeros(size(dists,2),1);
    for i=1:size(dists,2)
        ff=dists(:,i);
        pts=find(ff==min(ff));
        map(i)=start_points(pts);
    end

    assign=map;
    D=zeros(size(assign));
    uu=start_points;
    bds_total=[];
    for i=1:numel(uu)
        pts=find(assign==uu(i));
        bds=[];
        for p=1:numel(pts)
            gg=assign(NNs{pts(p)});
            if numel(unique(gg))>1
                bds(end+1)=pts(p);
            end
        end
        hh=intersect(pts,edges);
        bds=[bds(:);hh(:)];
        D(pts)=min(dist(pts,bds),[],2);
        bds_total(end+1:end+numel(bds))=bds;
    end
 end