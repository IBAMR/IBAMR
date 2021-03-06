%% ---------------------------------------------------------------------
%%
%% Copyright (c) 2014 - 2019 by the IBAMR developers
%% All rights reserved.
%%
%% This file is part of IBAMR.
%%
%% IBAMR is free software and is distributed under the 3-clause BSD
%% license. The full text of the license can be found in the file
%% COPYRIGHT at the top level directory of IBAMR.
%%
%% ---------------------------------------------------------------------

u32  = reshape(load('u.N=32' ), 33, 32);
u64  = reshape(load('u.N=64' ), 65, 64);
u128 = reshape(load('u.N=128'),129,128);
u256 = reshape(load('u.N=256'),257,256);

uI32  = 0.5* u64(1:2:size( u64,1),1:2:size( u64,2)) + 0.5* u64(1:2:size( u64,1),2:2:size( u64,2));
uI64  = 0.5*u128(1:2:size(u128,1),1:2:size(u128,2)) + 0.5*u128(1:2:size(u128,1),2:2:size(u128,2));
uI128 = 0.5*u256(1:2:size(u256,1),1:2:size(u256,2)) + 0.5*u256(1:2:size(u256,1),2:2:size(u256,2));

u1_32_64 = L1normx(u32-uI32);
u1_64_128 = L1normx(u64-uI64);
u1_128_256 = L1normx(u128-uI128);

fprintf('|u32 -u64 |_1  = %1.2e  rate = %2.2f\n',u1_32_64  ,log2(u1_32_64 /u1_64_128 ));
fprintf('|u64 -u128|_1  = %1.2e  rate = %2.2f\n',u1_64_128 ,log2(u1_64_128/u1_128_256));
fprintf('|u128-u256|_1  = %1.2e  rate = ---\n\n',u1_128_256);

u2_32_64 = L2normx(u32-uI32);
u2_64_128 = L2normx(u64-uI64);
u2_128_256 = L2normx(u128-uI128);

fprintf('|u32 -u64 |_2  = %1.2e  rate = %2.2f\n',u2_32_64  ,log2(u2_32_64 /u2_64_128 ));
fprintf('|u64 -u128|_2  = %1.2e  rate = %2.2f\n',u2_64_128 ,log2(u2_64_128/u2_128_256));
fprintf('|u128-u256|_2  = %1.2e  rate = ---\n\n',u2_128_256);

uoo_32_64 = max(max(abs(u32-uI32)));
uoo_64_128 = max(max(abs(u64-uI64)));
uoo_128_256 = max(max(abs(u128-uI128)));

fprintf('|u32 -u64 |_oo = %1.2e  rate = %2.2f\n',uoo_32_64  ,log2(uoo_32_64 /uoo_64_128 ));
fprintf('|u64 -u128|_oo = %1.2e  rate = %2.2f\n',uoo_64_128 ,log2(uoo_64_128/uoo_128_256));
fprintf('|u128-u256|_oo = %1.2e  rate = ---\n\n',uoo_128_256);

v32  = reshape(load('v.N=32' ), 32, 33);
v64  = reshape(load('v.N=64' ), 64, 65);
v128 = reshape(load('v.N=128'),128,129);
v256 = reshape(load('v.N=256'),256,257);

vI32  = 0.5* v64(1:2:size( v64,1),1:2:size( v64,2)) + 0.5* v64(2:2:size( v64,1),1:2:size( v64,2));
vI64  = 0.5*v128(1:2:size(v128,1),1:2:size(v128,2)) + 0.5*v128(2:2:size(v128,1),1:2:size(v128,2));
vI128 = 0.5*v256(1:2:size(v256,1),1:2:size(v256,2)) + 0.5*v256(2:2:size(v256,1),1:2:size(v256,2));

v1_32_64 = L1normy(v32-vI32);
v1_64_128 = L1normy(v64-vI64);
v1_128_256 = L1normy(v128-vI128);

fprintf('|v32 -v64 |_1  = %1.2e  rate = %2.2f\n',v1_32_64  ,log2(v1_32_64 /v1_64_128 ));
fprintf('|v64 -v128|_1  = %1.2e  rate = %2.2f\n',v1_64_128 ,log2(v1_64_128/v1_128_256));
fprintf('|v128-v256|_1  = %1.2e  rate = ---\n\n',v1_128_256);

v2_32_64 = L2normy(v32-vI32);
v2_64_128 = L2normy(v64-vI64);
v2_128_256 = L2normy(v128-vI128);

fprintf('|v32 -v64 |_2  = %1.2e  rate = %2.2f\n',v2_32_64  ,log2(v2_32_64 /v2_64_128 ));
fprintf('|v64 -v128|_2  = %1.2e  rate = %2.2f\n',v2_64_128 ,log2(v2_64_128/v2_128_256));
fprintf('|v128-v256|_2  = %1.2e  rate = ---\n\n',v2_128_256);

voo_32_64 = max(max(abs(v32-vI32)));
voo_64_128 = max(max(abs(v64-vI64)));
voo_128_256 = max(max(abs(v128-vI128)));

fprintf('|v32 -v64 |_oo = %1.2e  rate = %2.2f\n',voo_32_64  ,log2(voo_32_64 /voo_64_128 ));
fprintf('|v64 -v128|_oo = %1.2e  rate = %2.2f\n',voo_64_128 ,log2(voo_64_128/voo_128_256));
fprintf('|v128-v256|_oo = %1.2e  rate = ---\n\n',voo_128_256);

p32  = reshape(load('p.N=32' ), 32, 32);
p64  = reshape(load('p.N=64' ), 64, 64);
p128 = reshape(load('p.N=128'),128,128);
p256 = reshape(load('p.N=256'),256,256);

pI32  = 0.25* p64(1:2:size( p64,1),1:2:size( p64,1)) + 0.25* p64(2:2:size( p64,1),1:2:size( p64,1)) + 0.25* p64(1:2:size( p64,1),2:2:size( p64,1)) + 0.25* p64(2:2:size( p64,1),2:2:size( p64,1));
pI64  = 0.25*p128(1:2:size(p128,1),1:2:size(p128,1)) + 0.25*p128(2:2:size(p128,1),1:2:size(p128,1)) + 0.25*p128(1:2:size(p128,1),2:2:size(p128,1)) + 0.25*p128(2:2:size(p128,1),2:2:size(p128,1));
pI128 = 0.25*p256(1:2:size(p256,1),1:2:size(p256,1)) + 0.25*p256(2:2:size(p256,1),1:2:size(p256,1)) + 0.25*p256(1:2:size(p256,1),2:2:size(p256,1)) + 0.25*p256(2:2:size(p256,1),2:2:size(p256,1));

p1_32_64 = L1norm(p32-pI32);
p1_64_128 = L1norm(p64-pI64);
p1_128_256 = L1norm(p128-pI128);

fprintf('|p32 -p64 |_1  = %1.2e  rate = %2.2f\n',p1_32_64  ,log2(p1_32_64 /p1_64_128 ));
fprintf('|p64 -p128|_1  = %1.2e  rate = %2.2f\n',p1_64_128 ,log2(p1_64_128/p1_128_256));
fprintf('|p128-p256|_1  = %1.2e  rate = ---\n\n',p1_128_256);

p2_32_64 = L2norm(p32-pI32);
p2_64_128 = L2norm(p64-pI64);
p2_128_256 = L2norm(p128-pI128);

fprintf('|p32 -p64 |_2  = %1.2e  rate = %2.2f\n',p2_32_64  ,log2(p2_32_64 /p2_64_128 ));
fprintf('|p64 -p128|_2  = %1.2e  rate = %2.2f\n',p2_64_128 ,log2(p2_64_128/p2_128_256));
fprintf('|p128-p256|_2  = %1.2e  rate = ---\n\n',p2_128_256);

poo_32_64 = max(max(abs(p32-pI32)));
poo_64_128 = max(max(abs(p64-pI64)));
poo_128_256 = max(max(abs(p128-pI128)));

fprintf('|p32 -p64 |_oo = %1.2e  rate = %2.2f\n',poo_32_64  ,log2(poo_32_64 /poo_64_128 ));
fprintf('|p64 -p128|_oo = %1.2e  rate = %2.2f\n',poo_64_128 ,log2(poo_64_128/poo_128_256));
fprintf('|p128-p256|_oo = %1.2e  rate = ---\n\n',poo_128_256);
