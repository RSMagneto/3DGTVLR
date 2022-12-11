function h=fspecial_3D(type, f)

h=zeros(f,f,f);
switch type
    case 'disk'
        if mod(f,2)
            rad   =  floor(f / 2);
            crad  = ceil(rad-0.5);
            [x,y,z] = meshgrid(-crad:crad,-crad:crad,-crad:crad);
            xyz=sqrt(x.^2+y.^2+z.^2);
            ind= xyz<=rad;
            h(ind)=1;
        else
            rad   =  floor(f / 2);
            crad  = rad-0.5;
            [x,y,z] = meshgrid(-crad:crad,-crad:crad,-crad:crad);
            xyz=sqrt(x.^2+y.^2+z.^2);
            ind= xyz<=rad;
            h(ind)=1;
        end
        
end
end





