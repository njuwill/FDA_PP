
function [V_xb, V_yb,V_zb]=change_dir(grid_length,V_xb,V_yb,V_zb)
V_xb(:,end)=(2*double(V_xb(round(grid_length/2),end)>0)-1)*V_xb(:,end);
V_xb(:,end-1)=(2*double(V_xb(1,end-1)>0)-1)*V_xb(:,end-1);
V_yb(:,end)=(2*double(V_yb(round(grid_length/2),end)>0)-1)*V_yb(:,end);
V_yb(:,end-1)=(2*double(V_yb(round(grid_length/4),end-1)>0)-1)*V_yb(:,end-1);
V_zb(:,end)=(2*double(V_zb(round(grid_length/2),end)>0)-1)*V_zb(:,end);
V_zb(:,end-1)=(2*double(V_zb(round(grid_length/8),end-1)>0)-1)*V_zb(:,end-1);
