function pass = test_eigs_foxli(pref)
% Eigenvalues of the Fox-Li integral operator
% Toby Driscoll and Nick Trefethen, 7 October 2010

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-10;


% In the field of optics, integral operators arise that have a complex 
% symmetric (but not Hermitian) oscillatory kernel.  An example is the 
% following linear Fredholm operator L associated with the names of Fox 
% and Li (also Fresnel and H. J. Landau):
%
%    v(x) = sqrt(i*F/pi) int_{-1}^1 K(x,s) u(s) ds.
%
% L maps a function u defined on [-1,1] to another function v = Lu defined 
% on [-1,1].  The number F is a positive real parameter, the Fresnel number, 
% and the kernel function K(x,s) is 
%
%    K(x,s) = exp(-i*F*(x-s)^2).
%

% To create the operator in Chebfun, we define the kernel anduse the FRED
% function to build L:
F = 64*pi;                                     % Fresnel number
K = @(x,s) exp(-1i*F*(x-s).^2 );               % kernel
L = sqrt(1i*F/pi) * chebop(@(u) fred(K,u));    % Fredholm integral operator 

lam = eigs(L, 80, 'lm', pref);

lamV4 = [...
  0.999875850642719 + 0.002939727332351i
  0.999462827025035 + 0.011770704208759i
  0.998565291177526 + 0.026434202408201i
  0.997033328793154 + 0.046989672193378i
  0.994267953619479 + 0.073267872285422i
  0.990238812609936 + 0.105259125576187i
  0.983733882416718 + 0.142952367220762i
  0.974894475345593 + 0.185624643548422i
  0.962363797737049 + 0.234153430289937i
  0.945277602773947 + 0.286398628240589i
  0.924161315680008 + 0.343726584258937i
  0.895107269271094 + 0.404652622052363i
  0.861076446120097 + 0.466734538513570i
  0.818752223163400 + 0.534049646793061i
  0.765380225043449 + 0.597881327759312i
  0.708198659515296 + 0.662623141295791i
  0.634432079320151 + 0.728133518803939i
  0.466608483566232 + 0.836806328123111i
  0.554000245990309 + 0.781601712083035i
  0.360557208692383 + 0.880354558254032i
  0.252485785223710 + 0.909285660073678i
  0.135687620239551 + 0.931453254248935i
  0.008707591654773 + 0.935539084904055i
 -0.120665041568483 + 0.920788517860570i
 -0.250121395618950 + 0.885904471298700i
 -0.374957043710881 + 0.834229573916335i
 -0.501940083630190 + 0.762788028396196i
 -0.615406110737749 + 0.658336901566417i
 -0.704552083741678 + 0.542410547759310i
 -0.784288750050251 + 0.413116692989699i
 -0.842680161309719 + 0.260699287917244i
 -0.868976131554211 + 0.091972175632861i
 -0.851924611628668 - 0.077321651105772i
 -0.752018515704448 - 0.391798531165885i
 -0.812343247568437 - 0.230768551364885i
 -0.639557494409902 - 0.545682194967504i
 -0.491615488526506 - 0.660530053387882i
 -0.327389998513825 - 0.736111085640646i
 -0.156111451199873 - 0.779659477795052i
  0.030060156733991 - 0.791774914912482i
  0.224383335951275 - 0.751917075311678i
  0.401813223891073 - 0.659111115064116i
  0.550142763156226 - 0.521172166182755i
  0.651653259272269 - 0.343870505521980i
  0.694343859156978 - 0.154845355787791i
  0.692916870552612 + 0.030153892501502i
  0.650196338701473 + 0.212659062505210i
  0.556123330709173 + 0.389792447151947i
  0.392252454367271 + 0.530897352205454i
  0.198060039033987 + 0.587135127874017i
  0.027675611100829 + 0.583831636886728i
 -0.127961968960236 + 0.556672130288178i
 -0.290411438574085 + 0.490825381769936i
 -0.424576079907327 + 0.357861106045766i
 -0.496585678196401 + 0.199036572604063i
 -0.527654837806502 + 0.028594924690102i
 -0.492293349708395 - 0.159150442548895i
 -0.382028776132266 - 0.306101055610133i
 -0.241609868550439 - 0.397747976096999i
 -0.081755362113186 - 0.434933123173213i
  0.068768607599505 - 0.409597458822737i
  0.190133230461771 - 0.347253683518753i
  0.291717158538372 - 0.257784422268085i
  0.363438064569636 - 0.129186380134241i
  0.375848770921531 + 0.021881323583499i
  0.327773784870397 + 0.166667557801527i
  0.217245467984068 + 0.280044865624752i
  0.068023869101186 + 0.324636242377690i
 -0.073862016065679 + 0.288559972186171i
 -0.151432410914520 + 0.203137327813510i
 -0.232374973415348 + 0.061663343908047i
 -0.232584218812971 - 0.052805867120566i
 -0.190695410548372 + 0.136052058856451i
 -0.169954045481035 - 0.136970448685688i
 -0.095085668534019 - 0.181246288801342i
 -0.009779411172964 - 0.200011973310063i
  0.081020770717980 - 0.175426125911518i
  0.142262210069373 - 0.109679335395640i
  0.164615051419544 - 0.033865845184579i
  0.157090152579831 + 0.043277925117789i];


lam = sort(lam);
lamV4 = sort(lamV4);

err = norm(lam - lamV4, inf);
pass = err < tol;

end
