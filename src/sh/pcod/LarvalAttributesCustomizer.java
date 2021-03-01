/*
 * NonEggStageAttributesCustomizer.java
 *
 * Created on January 12, 2006, 4:20 PM
 */

package sh.pcod;

import wts.models.DisMELS.framework.*;
import wts.models.DisMELS.gui.AttributesCustomizer;

/**
 * @author William Stockhausen
 */
public class LarvalAttributesCustomizer extends AttributesCustomizer {

    private boolean showHorizPos = true;
    private boolean showVertPos = true;
    
    private AbstractLarvalAttributes attributes = null;
    
    /**
     * Creates new customizer NonEggStageAttributesCustomizer
     */
    public LarvalAttributesCustomizer() {
        initComponents();
    }
    
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jPanel1 = new javax.swing.JPanel();
        czrStandardAttributes = new wts.models.DisMELS.gui.AbstractLHSAttributesCustomizer();
        jPanel2 = new javax.swing.JPanel();
        jchkAttached = new javax.swing.JCheckBox();
        jtfDiameter = new javax.swing.JTextField();
        jLabel2 = new javax.swing.JLabel();

        setLayout(new java.awt.BorderLayout());

        jPanel1.setLayout(new java.awt.BorderLayout());
        jPanel1.add(czrStandardAttributes, java.awt.BorderLayout.CENTER);

        add(jPanel1, java.awt.BorderLayout.NORTH);

        jPanel2.setBorder(javax.swing.BorderFactory.createTitledBorder("Additional attributes"));

        jchkAttached.setText("attached?");
        jchkAttached.setActionCommand("");
        jchkAttached.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jchkAttachedItemStateChanged(evt);
            }
        });

        jtfDiameter.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
        jtfDiameter.setText("0");
        jtfDiameter.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jtfDiameterActionPerformed(evt);
            }
        });

        jLabel2.setText("length (mm)");
        jLabel2.setToolTipText("diameter in mm");

        org.jdesktop.layout.GroupLayout jPanel2Layout = new org.jdesktop.layout.GroupLayout(jPanel2);
        jPanel2.setLayout(jPanel2Layout);
        jPanel2Layout.setHorizontalGroup(
            jPanel2Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(jPanel2Layout.createSequentialGroup()
                .addContainerGap()
                .add(jchkAttached)
                .add(0, 0, Short.MAX_VALUE))
            .add(jPanel2Layout.createSequentialGroup()
                .add(jtfDiameter, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 123, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(jLabel2, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 259, Short.MAX_VALUE))
        );
        jPanel2Layout.setVerticalGroup(
            jPanel2Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(jPanel2Layout.createSequentialGroup()
                .add(jchkAttached)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(jPanel2Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jtfDiameter, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jLabel2))
                .addContainerGap())
        );

        add(jPanel2, java.awt.BorderLayout.CENTER);
    }// </editor-fold>//GEN-END:initComponents

    private void jtfDiameterActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jtfDiameterActionPerformed
        Double n = new Double(jtfDiameter.getText());
        attributes.setValue(attributes.PROP_SL,n);
    }//GEN-LAST:event_jtfDiameterActionPerformed

    private void jchkAttachedItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jchkAttachedItemStateChanged
        attributes.setValue(attributes.PROP_attached, jchkAttached.isSelected());
    }//GEN-LAST:event_jchkAttachedItemStateChanged

    @Override
    public void setObject(Object bean) {
        if (bean instanceof AbstractLarvalAttributes) {
            setAttributes((AbstractLarvalAttributes) bean);
        }
    }
    
    @Override
    public AbstractLarvalAttributes getAttributes() {
        return attributes;
    }
    
    @Override
    public void setAttributes(LifeStageAttributesInterface newAtts) {
        if (newAtts instanceof AbstractLarvalAttributes) {
            attributes = (AbstractLarvalAttributes) newAtts;
            czrStandardAttributes.setObject(attributes);
            Double d = null;
            jchkAttached.setSelected(attributes.getValue(attributes.PROP_attached,true));
            jtfDiameter.setText(attributes.getValue(attributes.PROP_SL,d).toString());
        }
    }
    
    @Override
    public void showHorizPos(boolean b) {
        czrStandardAttributes.showHorizPos(b);
        revalidate();
    }
    
    @Override
    public void showVertPos(boolean b) {
        czrStandardAttributes.showVertPos(b);
        revalidate();
    }
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private wts.models.DisMELS.gui.AbstractLHSAttributesCustomizer czrStandardAttributes;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JCheckBox jchkAttached;
    private javax.swing.JTextField jtfDiameter;
    // End of variables declaration//GEN-END:variables
    
}