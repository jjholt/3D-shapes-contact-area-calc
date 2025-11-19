function d = pointTriangleDistance(p, a, b, c)
% pointTriangleDistance
% ---------------------
% Exact shortest distance from point p to triangle (a,b,c).
% Based on standard closest-point-on-triangle algorithm.

% Edges
ab = b - a;
ac = c - a;
ap = p - a;

d1 = dot(ab, ap);
d2 = dot(ac, ap);

% Check if P in vertex region outside A
if d1 <= 0 && d2 <= 0
    d = norm(p - a);
    return;
end

% Check if P in vertex region outside B
bp = p - b;
d3 = dot(ab, bp);
d4 = dot(ac, bp);
if d3 >= 0 && d4 <= d3
    d = norm(p - b);
    return;
end

% Check if P in edge region of AB
vc = d1*d4 - d3*d2;
if vc <= 0 && d1 >= 0 && d3 <= 0
    v = d1 / (d1 - d3);
    proj = a + v * ab;
    d = norm(p - proj);
    return;
end

% Check if P in vertex region outside C
cp = p - c;
d5 = dot(ab, cp);
d6 = dot(ac, cp);
if d6 >= 0 && d5 <= d6
    d = norm(p - c);
    return;
end

% Check if P in edge region of AC
vb = d5*d2 - d1*d6;
if vb <= 0 && d2 >= 0 && d6 <= 0
    w = d2 / (d2 - d6);
    proj = a + w * ac;
    d = norm(p - proj);
    return;
end

% Check if P in edge region of BC
va = d3*d6 - d5*d4;
if va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0
    w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    proj = b + w * (c - b);
    d = norm(p - proj);
    return;
end

% P inside face region
denom = 1 / (va + vb + vc);
v = vb * denom;
w = vc * denom;
proj = a + ab * v + ac * w;
d = norm(p - proj);

end