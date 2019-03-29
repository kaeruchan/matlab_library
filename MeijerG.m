function out = MeijerG( an, ap, bm, bq, z )

an_str = vec2str(an);
ap_str = vec2str(ap);
bm_str = vec2str(bm);
bq_str = vec2str(bq);
z_str = num2str(z, 32);

MeijerGmupad_str = ['float(meijerG(', an_str, ',', ap_str, ',', ...
    bm_str, ',', bq_str, ',' z_str, '))'];

out = double(evalin(symengine, MeijerGmupad_str)); 

return;


function an_str = vec2str(an)

if isempty(an)
    an_str = '[]';
else
    an_str = ['[', num2str(an(1), 32)];
    for i=2:length(an)
        an_str = [an_str, ', ', num2str(an(i), 32)];
    end;
    an_str = [an_str, ']'];
end;

return;