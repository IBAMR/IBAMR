# Visualizing IBAMR output

## With VisIt

As of IBAMR 0.19, the visualization output files no longer contain separate
`U_x`, `U_y`, and `U_z` variables: instead, all velocity values are represented
by the vector-valued `U` field. To visualize independent components, you can add
new variables in the Expressions editor. For example, to get the first
component, define a new variable `Ux` as `U[0]`.

We recommend setting up a consistent environment for VisIt variables so that all
quantities are always defined. VisIt saves variables in xml files: here is one
such file which defines all velocity components:

```xml
<?xml version="1.0"?>
<Object name="ExpressionList">
    <Object name="Expression">
        <Field name="name" type="string">Ux</Field>
        <Field name="definition" type="string">U[0]</Field>
        <Field name="hidden" type="bool">false</Field>
        <Field name="type" type="string">ScalarMeshVar</Field>
        <Field name="fromDB" type="bool">false</Field>
        <Field name="fromOperator" type="bool">false</Field>
        <Field name="operatorName" type="string">__none__</Field>
        <Field name="meshName" type="string"></Field>
        <Field name="dbName" type="string">__none__</Field>
        <Field name="autoExpression" type="bool">false</Field>
    </Object>
    <Object name="Expression">
        <Field name="name" type="string">Uy</Field>
        <Field name="definition" type="string">U[1]</Field>
        <Field name="hidden" type="bool">false</Field>
        <Field name="type" type="string">ScalarMeshVar</Field>
        <Field name="fromDB" type="bool">false</Field>
        <Field name="fromOperator" type="bool">false</Field>
        <Field name="operatorName" type="string">__none__</Field>
        <Field name="meshName" type="string"></Field>
        <Field name="dbName" type="string">__none__</Field>
        <Field name="autoExpression" type="bool">false</Field>
    </Object>
    <Object name="Expression">
        <Field name="name" type="string">Uz</Field>
        <Field name="definition" type="string">U[2]</Field>
        <Field name="hidden" type="bool">false</Field>
        <Field name="type" type="string">ScalarMeshVar</Field>
        <Field name="fromDB" type="bool">false</Field>
        <Field name="fromOperator" type="bool">false</Field>
        <Field name="operatorName" type="string">__none__</Field>
        <Field name="meshName" type="string"></Field>
        <Field name="dbName" type="string">__none__</Field>
        <Field name="autoExpression" type="bool">false</Field>
    </Object>
</Object>
```
