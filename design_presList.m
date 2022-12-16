function [pList] = design_presList(dt,lon,lat,varargin)

st = -20;
p_tol = 1;
tol_m = 1.01;
scena = 1;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'STEP'
                st = varargin{i+1};
            case 'P_TOL'
                p_tol = varargin{i+1};
            case 'TOL_M'
                tol_m = varargin{i+1};
            case 'SCENA'
                scena = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

[param0] = mcd_get_vertical_profile(dt,lon,lat,'XZ',0,'verbose',false,'scena',scena);

pList = [(param0.pres-1e-4)];
while(pList(end)>p_tol)
    pList_tmp = pList(end)+st;
    pList_tmp2 = pList(end)+st*2;
    
    param_end = mcd_get_vertical_profile(dt,lon,lat,'XZ',pList(end),'zkey',4,'verbose',false,'scena',scena);
    [param_tmp] = mcd_get_vertical_profile(dt,lon,lat,'XZ',pList_tmp,'zkey',4,'verbose',false,'scena',scena);
    [param_tmp2] = mcd_get_vertical_profile(dt,lon,lat,'XZ',pList_tmp2,'zkey',4,'verbose',false,'scena',scena);
    
    
    if param_tmp.alt_above_surf*tol_m < (param_end.alt_above_surf+param_tmp2.alt_above_surf)/2
        st = st / 2;
    else
        pList = [pList pList_tmp];
    end
end