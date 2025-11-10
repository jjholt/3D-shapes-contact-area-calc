---

## Bounding Box Distance Method

The updated contact detection method treats the cube as a **solid, axis-aligned volume**, rather than a collection of surface vertices.
For each triangle on the pyramid mesh, contact is determined by measuring how close its centroid lies to the cube’s bounding box.

---

### **Step 1 — Compute Triangle Centroid and Area**

Each face has three vertices.
The centroid is simply the average of those vertices:

```math
C = \frac{p_1 + p_2 + p_3}{3}
```

---

### **Step 2 — Define the Cube’s Bounding Box**

From all cube vertices `Vc`, determine the minimum and maximum coordinates:

```text
xmin, xmax
ymin, ymax
zmin, zmax
```

These define the smallest rectangular box that encloses the cube.
Any point `(x, y, z)` inside this box satisfies:

```text
xmin ≤ x ≤ xmax
ymin ≤ y ≤ ymax
zmin ≤ z ≤ zmax
```

---

### **Step 3 — Find the Closest Point on the Box to the Centroid**

For each centroid coordinate `(cx, cy, cz)`, **clamp** it to the cube bounds:

```text
qx = clamp(cx, xmin, xmax)
qy = clamp(cy, ymin, ymax)
qz = clamp(cz, zmin, zmax)
```

* If the centroid is **inside** the cube along that axis → the coordinate is unchanged.
* If it’s **outside**, it’s clamped to the nearest cube face.

The resulting point
`Q = (qx, qy, qz)`
is the **closest point on or inside the cube** to the centroid.

---

### **Step 4 — Compute Distance to the Closest Point**

```math
d = \sqrt{(cx - qx)^2 + (cy - qy)^2 + (cz - qz)^2}
```

* Centroid inside cube → `d = 0`
* Centroid just outside cube → `d` equals perpendicular distance to the nearest cube face
* Centroid far from cube → `d` increases with separation

---

### **Step 5 — Mark Contact if Within Tolerance**

A triangle is considered “in contact” if:

```text
d ≤ tol
```

This means:

* ✅ **Embedded or touching triangles** → `d ≈ 0`
* ✅ **Near-contact triangles** → within tolerance band
* ❌ **Far-away triangles** → ignored

---

This approach robustly identifies **flat-on-flat** contact areas (like a pyramid base resting on a cube), avoids missing faces due to sparse vertex proximity, and allows simple control of sensitivity through the tolerance parameter `tol`.
