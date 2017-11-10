/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package OMEuiUtils;



import Glacier2.CannotCreateSessionException;
import Glacier2.PermissionDeniedException;
import java.awt.Dimension;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.BorderFactory;
import javax.swing.JButton;


import java.awt.BorderLayout;
import java.awt.Dialog;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;
import java.util.Map;
import java.util.Scanner;
import javax.swing.DefaultListModel;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.JTree;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
//import javax.swing.event.TreeExpansionEvent;
//import javax.swing.event.TreeWillExpandListener;

import javax.swing.tree.DefaultMutableTreeNode;
//import javax.swing.tree.ExpandVetoException;

import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;

import omero.ServerError;
import omero.api.IContainerPrx;
import omero.api.IMetadataPrx;
import omero.api.ServiceFactoryPrx;

import omero.client;
import omero.model.Dataset;
import omero.model.Experimenter;
import omero.model.IObject;
import omero.model.Image;
import omero.model.Plate;
import omero.model.Project;
import omero.model.Screen;
import omero.sys.ParametersI;
import omero.gateway.model.DatasetData;
import omero.gateway.model.ImageData;
import omero.gateway.model.PlateData;
import omero.gateway.model.ProjectData;
import omero.gateway.model.ScreenData;
import omero.model.FileAnnotation;
import omero.model.OriginalFile;



/**
 *
 * @author imunro
 */



public class OMEROImageChooser extends JDialog implements ActionListener {

    
    private JTree tree;
    private ArrayList<Object> returned; 
    private TreePath expPath;
    private JTextField fnameField; 
    private JComboBox extList;
    private boolean promptForInput;
    // 0 == image, 1== dataset, 2 = Plate, 
    // NB Not selectable   3 == project, 4 = Screen, 5 = user
    // 6 = File Attachment to Dataset, 7 = File Attachment to Plate
    // 8 = Image OR File Attachment to Dataset
    private final int selectedType;
    
  
    // Used to get Attachments
    private JList listbox;
    private DefaultListModel attachmentNameModel;
    private ArrayList<FileAnnotation> attachmentList;
    private ArrayList<String> annotationType;
    private IMetadataPrx metadataService;
    private ParametersI attachmentParam;
    private OriginalFile returnedFile;
    
    
    private boolean loadImages = false;
    
    
    // Select dataset for output
    public OMEROImageChooser(omero.client omeroclient, long userId, Long expandId, String[] filenameStrings)  {
      this(omeroclient, userId, 1, false, expandId, filenameStrings);
    }
    
    // Select any object by type
    public OMEROImageChooser(omero.client omeroclient, long userId ,int selectedType)  {
      this(omeroclient, userId, selectedType, false, new Long(-1), null);
    }
    
    // allow selection of multiple Images
    public OMEROImageChooser(omero.client omeroclient, long userId, boolean allowMultiple )  {
      this(omeroclient, userId, 0, allowMultiple, new Long(-1), null);
    }
    
    // select a single image 
    public OMEROImageChooser(omero.client omeroclient,  long userId, Long expandId)  {
      this(omeroclient, userId, 0, false, expandId, null);
    }
    
    // Expand dataset (type = 0 single Image , 6 attachment to Dataset or 8 Image-or-Attachment 
    public OMEROImageChooser(omero.client omeroclient,  long userId, int type, Long expandId)  {
      this(omeroclient, userId, type, false, expandId, null);
    }
    
    // Expand dataset (Images) 
    public OMEROImageChooser(omero.client omeroclient, long userId, boolean allowMultiple, Long expandId)  {
      this(omeroclient, userId, 0,  allowMultiple,  expandId, null);
    }
    

    public OMEROImageChooser(omero.client omeroclient, long userId,  int type, boolean allowMultiple,  Long expandId, String[] filenameStrings)  {
           
      
      returned = null;
      returnedFile = null;
      
      expPath = null;
     
      setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
      
      // Create a Dialog to display while tree is populated
      final JDialog waitDialog = new JDialog(this, "Fetching data. Please Wait...", false);
    
      waitDialog.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
      //waitDialog.setUndecorated(true);
      waitDialog.pack();
      waitDialog.setSize(400,60);
      waitDialog.setLocationRelativeTo(null);
  
      waitDialog.setVisible(true); // show the dialog on the screen
      
      
      String name;
      try {
        
        this.setModalityType(Dialog.ModalityType.APPLICATION_MODAL);
        //ServiceFactoryPrx session = omeroclient.joinSession(sessionid);
        ServiceFactoryPrx session = omeroclient.getSession();
        
        IContainerPrx proxy = session.getContainerService();
        ParametersI param = new ParametersI();
        ParametersI paramAll = new ParametersI();
       
        param.exp(omero.rtypes.rlong(userId));
        paramAll.exp(omero.rtypes.rlong(userId));
       
        System.out.println("Uid = " + userId);
        Experimenter exp = session.getAdminService().getExperimenter(userId);
        name = exp.getFirstName().getValue() + " " + exp.getLastName().getValue();
        
       // name = session.getAdminService().getEventContext().userName;
        datasetInfo userInfo = new datasetInfo(name, 0L, 5); // type 5 is a user Id 0
        DefaultMutableTreeNode userNode = new DefaultMutableTreeNode(userInfo);

        //create the tree by passing in the root node
        tree = new JTree(userNode);
        tree.setShowsRootHandles(true);
        //tree.setRootVisible( false );
        
                    
        JScrollPane spane = new JScrollPane(tree);
        spane.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), " ") );
        spane.setMinimumSize(new Dimension(350,300));
        spane.setPreferredSize(new Dimension(350,600));
        
        spane.add(tree);
        spane.setViewportView(tree);
        
        JPanel buttonPanel = new JPanel();
        
        final JButton openButton = new JButton("Open");
        openButton.setActionCommand("Open");
        openButton.setEnabled(false);
        
        attachmentNameModel = new DefaultListModel();
        attachmentList = new ArrayList<>();
        
        boolean isPlate = false;  //default
        
        switch (type) {
          case 1:
            tree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
            if (promptForInput) {
              setTitle("Please select a Dataset, " + filenameStrings[1] + " & type");
            } else {
              setTitle("Please select a Dataset");
            }
            //param.noLeaves(); //no images loaded, this is the default value.
            break;
          case 2:
            isPlate = true;
            tree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
            if (promptForInput) {
              setTitle("Please select a Plate, " + filenameStrings[1] + " & type");
            } else {
              setTitle("Please select a Plate");
            }
            //datasetsList = proj.linkedPlateList;
            //param.noLeaves(); //no images loaded, this is the default value.
            break;
          case 6:  // 6 = File Attachment to Dataset
            createAttachmentPanel(openButton, userId, session);
            tree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
            setTitle("Please select an Attachment");
             //param.noLeaves(); //no images loaded, this is the default value.
            break;
            
          case 7:   // 7 = File Attachment to Plate
            type = 6;
            createAttachmentPanel(openButton, userId, session);
            isPlate = true;
            tree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
            setTitle("Please select an Attachment");
            break;

          // type 8 is a hybrid. Allows the user to choose either an image or a File Annotation to a Dataset
          case 8:
            type = 6;
            createAttachmentPanel(openButton, userId, session);
            tree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
            setTitle("Please select an Image or Attachment");
            paramAll.leaves();  //indicate to load the images
            loadImages = true;
            break;

          default:
            if (allowMultiple) {
              tree.getSelectionModel().setSelectionMode(TreeSelectionModel.DISCONTIGUOUS_TREE_SELECTION);
              setTitle("Please select one or more Images");
            } else {
              tree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
              setTitle("Please select an Image");
            }
            // Add pre-expansion event listener
            //tree.addTreeWillExpandListener(new MyTreeWillExpandListener());
            loadImages = true;
            paramAll.leaves();  //indicate to load the images
            break;
        }
        
        
        //Listen for when the tree selection changes.
        tree.addTreeSelectionListener(new TreeSelectionListener() {
          @Override
          public void valueChanged(TreeSelectionEvent e) {
            DefaultMutableTreeNode node = (DefaultMutableTreeNode) e.getPath().getLastPathComponent();
            attachmentList.clear();
            attachmentNameModel.clear();
            openButton.setEnabled(false);  // disable by default
            
              if ( selectedType==6) {  // 6 ==FileAttachment 
                // Allow selection of an image rather than an Attachment where required
                if (node.isLeaf() & loadImages) {
                  openButton.setEnabled(true);
                }
                datasetInfo di = (datasetInfo) node.getUserObject();
                Long objId = new Long(-1);
                String parentType = null;
                // show dataset and plate attachments only ATM
                if (di.getType() == 1 ) {   
                  Dataset dataset = ((DatasetData) di.getObject()).asDataset();
                  objId = dataset.getId().getValue();
                  parentType = "omero.model.Dataset";
                }
                if (di.getType() == 2 ) {   
                  Plate plate = ((PlateData) di.getObject()).asPlate();
                  objId = plate.getId().getValue();
                  parentType = "omero.model.Plate";
                }
                if (objId != -1) {
                  ArrayList<Long> Ids = new ArrayList<>();
                  Ids.add(objId);
                  List<Long> annotators = null;
                  Map<Long, List<IObject>> map = null;

                  try {
                    map = metadataService.loadAnnotations(parentType, Ids, annotationType, annotators, attachmentParam);
                  } catch (ServerError ex) {
                    Logger.getLogger(OMEROImageChooser.class.getName()).log(Level.SEVERE, null, ex);
                  }
                  
                  List<IObject> annotations = map.get(objId);
                  for (int a = 0; a < annotations.size(); a++) {
                    IObject obj = annotations.get(a);
                    if (obj instanceof FileAnnotation) {
                      FileAnnotation ann = (FileAnnotation) obj;
                      attachmentList.add(ann);
                      attachmentNameModel.addElement(ann.getFile().getName().getValue());
                    }
                  }
                }
              }
            else // For all types other than 6 
            // Enable open button if a leaf selected
            if (node.isLeaf()) {
              if (!promptForInput) {
                openButton.setEnabled(true);
              } else if (!fnameField.getText().isEmpty()) {
                openButton.setEnabled(true);
              }
            }
          }
        });
       
       
        // if choosing a filename/s for output
        if (filenameStrings != null)  {
          promptForInput = true;
         
          fnameField = new JTextField(filenameStrings[0], 8);
          fnameField.getDocument().addDocumentListener(new DocumentListener() {
           
            public void check() {
              if (fnameField.getText().isEmpty()) {
                openButton.setEnabled(false);
              } else {
                openButton.setEnabled(true);
              }
            }

            
            @Override
            public void insertUpdate(DocumentEvent e) {
              if (!fnameField.getText().isEmpty()) {
                TreePath[] paths = tree.getSelectionPaths();
                ArrayList<Object> selected = new ArrayList<Object>();
      
                if (paths != null) {
                  for (TreePath path : paths) {
                    DefaultMutableTreeNode node = (DefaultMutableTreeNode) (path.getLastPathComponent());
                    if (node.isLeaf()) {
                      openButton.setEnabled(true);
                    }
                  }
                }
              }
            }

            @Override
            public void removeUpdate(DocumentEvent e) {
                if (fnameField.getText().isEmpty()) {
                  openButton.setEnabled(false);
                }
            }

            @Override
            public void changedUpdate(DocumentEvent e) {
              
            } 

          });
          
          buttonPanel.add(fnameField);

          String[] extStrings = Arrays.copyOfRange(filenameStrings, 2, filenameStrings.length);
          extList = new JComboBox(extStrings);
          buttonPanel.add(extList);
          openButton.setText("Save");
          
        }
        else {
          promptForInput = false;
        }
        
        
        JButton cancelButton = new JButton("Cancel");
        cancelButton.setActionCommand("Cancel");
      
        buttonPanel.add(cancelButton, BorderLayout.LINE_START);
        // register the ButtonFrame object as the listener for the JButton.
        cancelButton.addActionListener( this ); 
          
        buttonPanel.add(openButton, BorderLayout.LINE_END);
        // register the ButtonFrame object as the listener for the JButton.
        openButton.addActionListener( this ); 
       
        add(buttonPanel, BorderLayout.SOUTH );
        add(spane);
   
       
        
        
        if (isPlate)  {   // Plate requested
          
          
          List<IObject> screenList = proxy.loadContainerHierarchy(Screen.class.getName(), new ArrayList<Long>(), param);
          List<IObject> allplatesList = proxy.loadContainerHierarchy(Plate.class.getName(), new ArrayList<Long>(), param);

          Collections.sort(screenList, new Comparator<IObject>() {
            @Override
            public int compare(IObject projOne, IObject projTwo) {
              return (new ScreenData((Screen)projOne).getName().compareToIgnoreCase(new ScreenData((Screen)projTwo).getName()));
            }
          }); 
          
          
          Iterator<IObject> i = screenList.iterator();
          ScreenData screen;
          java.util.Set<PlateData> plates;
          Iterator<PlateData> j;
          PlateData plate;
          PlateData pl;
          while (i.hasNext()) {
            screen = new ScreenData((Screen) i.next());

            plates = screen.getPlates();

            String screenName = screen.getName() + " [" + Integer.toString(plates.size()) + "]";
            long sId = screen.getId();
            datasetInfo screenInfo = new datasetInfo(screenName, sId, 4);    //type 4 is a screen
            DefaultMutableTreeNode projNode = new DefaultMutableTreeNode(screenInfo);
            userNode.add(projNode);

            List<PlateData> plateList = new ArrayList<PlateData>(plates);
            Collections.sort(plateList, new Comparator<PlateData>() {
              @Override
              public int compare(PlateData pOne, PlateData pTwo) {
                return ( ((PlateData)pOne).getName().compareToIgnoreCase( ((PlateData)pTwo).getName()));
              }
            }); 

            j = plateList.iterator();
            while (j.hasNext()) {
              plate = j.next();
              long pId = plate.getId();
              for( int ap=0; ap < allplatesList.size(); ap++)  {
                if (allplatesList.get(ap).getId().getValue() == pId)  {
                  pl = new PlateData((Plate)allplatesList.get(ap));
                  DefaultMutableTreeNode node = addPlate(plate);
                  projNode.add(node);
                  if (pId == expandId) {
                    expPath = new TreePath(node.getPath());
                  }
                  allplatesList.remove(ap);
                  break;
                }
              }
            }  
            
          }
          Collections.sort(allplatesList, new Comparator<IObject>() {
            @Override
            public int compare(IObject plOne, IObject plTwo) {
              return (new PlateData((Plate)plOne).getName().compareToIgnoreCase(new PlateData((Plate)plTwo).getName()));
            }
          }); 
          

          for (int p = 0; p < allplatesList.size(); p++)  {
            pl = new PlateData((Plate)allplatesList.get(p));
            DefaultMutableTreeNode node =  addPlate(pl);
            userNode.add(node);  
          }
          
          

        }
        else {  
          
          List<IObject> projectList = proxy.loadContainerHierarchy(Project.class.getName(), new ArrayList<Long>(), param);
          List<IObject> alldatasetsList = proxy.loadContainerHierarchy(Dataset.class.getName(), new ArrayList<Long>(), paramAll);

          Collections.sort(projectList, new Comparator<IObject>() {
            @Override
            public int compare(IObject projOne, IObject projTwo) {
              return (new ProjectData((Project)projOne).getName().compareToIgnoreCase(new ProjectData((Project)projTwo).getName()));
            }
          }); 
          
          Iterator<IObject> i = projectList.iterator();
          ProjectData project;
          java.util.Set<DatasetData> datasets;
          Iterator<DatasetData> j;
          Iterator<IObject> k = alldatasetsList.iterator();
          DatasetData dataset;
          DatasetData dset;
          String projName;
          String projNamePlus;
          
          while (i.hasNext()) {
            project = new ProjectData((Project) i.next());

            datasets = project.getDatasets();
            projName = project.getName();
            projNamePlus = projName + " [" + Integer.toString(datasets.size()) + "]";
            long pId = project.getId();
            datasetInfo projInfo = new datasetInfo(projNamePlus, pId, 3);    //type 3 is a project
            DefaultMutableTreeNode projNode = new DefaultMutableTreeNode(projInfo);
            userNode.add(projNode);
            List<DatasetData> datasetList = new ArrayList<DatasetData>(datasets);
            Collections.sort(datasetList, new Comparator<DatasetData>() {
              @Override
              public int compare(DatasetData dOne, DatasetData dTwo) {
                return ( dOne.getName().compareToIgnoreCase( dTwo.getName()));
              }
            }); 

            j = datasetList.iterator();
            while (j.hasNext()) {
              dataset = j.next();
              long dId = dataset.getId();
              for( int ad=0; ad < alldatasetsList.size(); ad++)  {
                if (alldatasetsList.get(ad).getId().getValue() == dId)  {
                  dset = new DatasetData((Dataset)alldatasetsList.get(ad));
                  DefaultMutableTreeNode node = addDataset(dset);
                  projNode.add(node);
                  if (dId == expandId) {
                    expPath = new TreePath(node.getPath());
                  }
                  alldatasetsList.remove(ad);
                  break;
                }
              }
            }   
          }
          
          Collections.sort(alldatasetsList, new Comparator<IObject>() {
            @Override
            public int compare(IObject dsetOne, IObject dsetTwo) {
              return (new DatasetData((Dataset)dsetOne).getName().compareToIgnoreCase(new DatasetData((Dataset)dsetTwo).getName()));
            }
          }); 
          

          for (int d = 0; d < alldatasetsList.size(); d++)  {
            dset = new DatasetData((Dataset)alldatasetsList.get(d));
            DefaultMutableTreeNode node =  addDataset(dset);
            userNode.add(node);  
          }
        }

        customCellRenderer renderer = new customCellRenderer();
        renderer.setIcons();
        tree.setCellRenderer(renderer);

      } catch (ServerError ex) {
        Logger.getLogger(OMEROImageChooser.class.getName()).log(Level.SEVERE, null, ex);
      }

      this.selectedType = type;

      if (expPath != null) {
        tree.expandPath(expPath);
        tree.scrollPathToVisible(expPath);
        tree.setSelectionPath(expPath);
      } else {
        tree.expandRow(0);
      }

     
      waitDialog.setVisible(false);
      waitDialog.dispose();
      
      setSize(400,400);
      setLocationRelativeTo(null);
     
      pack();
      
      setVisible(true);
      
    }
    
    public Image[] getSelectedImages()  {
      
      if (loadImages & returned != null)  {
        
        return returned.toArray(new Image[returned.size()]);
      }
      else {
        return new Image[0];
      }
    }
    
    public Dataset getSelectedDataset()  {
      
      if (selectedType == 1 & returned != null)  {  
        return (Dataset)returned.get(0);  
      }
      else {
        return null;
      }   
    }
    
    public String getFilename()  {
      
      if ( returned != null & promptForInput)  {  
        String fname = fnameField.getText() + extList.getSelectedItem().toString();
        return fname;    
      }
      else {
        return null;
      }   
    }
    
    public Plate getSelectedPlate()  {
      
      if (selectedType == 2 & returned != null)  {   
        return (Plate)returned.get(0);
      }
      else {
        return null;
      }   
    }
    
    public OriginalFile getSelectedFile()  {
      
      if (selectedType == 6 & returnedFile != null)  {   
        return returnedFile;
      }
      else {
        return null;
      }   
    }
    
    
    private void createAttachmentPanel(final JButton openButton, long userId, ServiceFactoryPrx session) {
      
    ListSelectionListener attachmentListener = new ListSelectionListener() {
      @Override
      public void valueChanged(ListSelectionEvent e) {
        int selectedAttachment = listbox.getSelectedIndex();
        if (selectedAttachment > -1 & selectedAttachment < attachmentNameModel.size()) {
          openButton.setEnabled(true);
        }
      }
    };

    // Create a new listbox control
    listbox = new JList(attachmentNameModel);
    listbox.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    listbox.addListSelectionListener(attachmentListener);

    JScrollPane attachmentPane = new JScrollPane();
    attachmentPane.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "Attachments"));
    attachmentPane.setMinimumSize(new Dimension(200, 300));
    attachmentPane.setPreferredSize(new Dimension(200, 300));
    attachmentPane.add(listbox);
    attachmentPane.setViewportView(listbox);
    add(attachmentPane, BorderLayout.EAST);

    annotationType = new ArrayList<>();
    annotationType.add("ome.model.annotations.FileAnnotation");
    attachmentParam = new ParametersI();
    attachmentParam.exp(omero.rtypes.rlong(userId)); //load the annotation for a given user.
    try {
      metadataService = session.getMetadataService();
    } catch (ServerError ex) {
      Logger.getLogger(OMEROImageChooser.class.getName()).log(Level.SEVERE, null, ex);
    }
  }
 
    private DefaultMutableTreeNode addPlate(PlateData plate)    {

      String pName = plate.getName();
      long pId = plate.getId();
      datasetInfo dsetInfo = new datasetInfo(pName, pId, 2 );  // type 2 is a plate
      dsetInfo.setObject(plate);
      DefaultMutableTreeNode node  = new DefaultMutableTreeNode(dsetInfo);
      
      return node;
    }

    private DefaultMutableTreeNode addDataset(DatasetData dataset)    {

      
      java.util.Set<ImageData> images = null;
      
      String dsetName = dataset.getName();
      if (loadImages)  {
         images = dataset.getImages();
        dsetName += " [" + Integer.toString(images.size()) + "]";
      }

      long dId = dataset.getId();
      datasetInfo dsetInfo = new datasetInfo(dsetName, dId, 1 );  // type 1 is a dataset
      dsetInfo.setObject(dataset);
      DefaultMutableTreeNode node  = new DefaultMutableTreeNode(dsetInfo);
      if (loadImages)  {
        node = addImages( dataset, node, images);
      }
      return node;
    }
    
    private DefaultMutableTreeNode addImages(DatasetData dataset, DefaultMutableTreeNode node, java.util.Set<ImageData> images)    {
      

      List<ImageData> imageList = new ArrayList<>(images);
      Collections.sort(imageList, new Comparator<ImageData>() {
        @Override
        public int compare(ImageData iOne, ImageData iTwo) {
          return (iOne.getName().compareToIgnoreCase(iTwo.getName()));
        }
      });
      
      Iterator<ImageData> j = imageList.iterator();
      ImageData image;
      while (j.hasNext()) {
        image = j.next();
        String imName = image.getName();
        Long imId = image.getId();
        datasetInfo imageInfo = new datasetInfo(imName, imId, 0 );  // type 0 is an image
        imageInfo.setObject(image);
        DefaultMutableTreeNode imNode = new DefaultMutableTreeNode(imageInfo);
        node.add(imNode);
      }
      return node;
      
    }

  @Override
  public void actionPerformed(ActionEvent e) {
    datasetInfo info = null;
    String command = e.getActionCommand();
    if (command.equals("Open")) {
      TreePath[] paths = tree.getSelectionPaths();
      ArrayList<Object> selected = new ArrayList<Object>();
      
      if (paths != null) {
        for (TreePath path : paths) {
          DefaultMutableTreeNode node = (DefaultMutableTreeNode) (path.getLastPathComponent());
          if (node.isLeaf() | selectedType == 6) {
            datasetInfo di = (datasetInfo) node.getUserObject();
            switch (selectedType) {
              case 1:
                if (di.getType() == 1) {
                  selected.add(((DatasetData) di.getObject()).asDataset());
                }
                break;
              case 2:
                if (di.getType() == 2) {
                  selected.add(((PlateData) di.getObject()).asPlate());
                }
                break;
              case 6:
                Integer ditype = di.getType();
                if (ditype == 1 | ditype == 2) {  
                  int selectedAttachment = listbox.getSelectedIndex();
                  if (selectedAttachment != -1) {   
                    FileAnnotation attachment = attachmentList.get(selectedAttachment);
                    String na = attachment.getFile().getName().getValue();
                    // Check names match just in case. Should always match.
                    if(na.equalsIgnoreCase((String)attachmentNameModel.getElementAt(selectedAttachment)))  {
                       returnedFile = attachment.getFile();
                    }
                  }
                }
                // also allow an image to be selected where appropriate 
                if (loadImages) {
                  if (di.getType() == 0) {
                    selected.add(((ImageData) di.getObject()).asImage());
                  }
                }
                break;
              default:
                if (di.getType() == selectedType) {
                  selected.add(((ImageData) di.getObject()).asImage());
                }
                break;
            }
          }

        }
      }
      if (!selected.isEmpty()) {
        returned = selected;
      }

      setVisible(false);
      dispose();

    }
    
    if (command.equals("Cancel")) {
      setVisible(false);
      dispose(); 
    }
  }
  
  
      
// Pre-expansion/collapse event listener
/*public class MyTreeWillExpandListener implements TreeWillExpandListener {
    public void treeWillExpand(TreeExpansionEvent evt) throws ExpandVetoException {
      JTree tree = (JTree) evt.getSource();

      // Get the path that will be expanded
      TreePath path = evt.getPath();
      DefaultMutableTreeNode node = (DefaultMutableTreeNode) path.getLastPathComponent();
      datasetInfo dInfo = (datasetInfo) node.getUserObject();
      int type = dInfo.getType();
      DefaultMutableTreeNode leaf = node.getFirstLeaf();
      if (type == 3) {     //Project
        for (int c = 0; c < node.getLeafCount(); c++) {
          datasetInfo lInfo = (datasetInfo) leaf.getUserObject();
          if (lInfo.getType() == 1) {  //Dataset
            DatasetData d = (DatasetData) lInfo.getObject();
            java.util.Set<ImageData> images = d.getImages();
            addImages(d, leaf, images);
          }
          leaf = node.getNextLeaf();
        } 
      } 


    } 

   public void treeWillCollapse(TreeExpansionEvent evt) throws ExpandVetoException {
        JTree tree = (JTree)evt.getSource();

        // Get the path that will be collapsed
        //TreePath path = evt.getPath();
        //Object e = evt.getSource();
       

        
    } 
}  */

    
    
    public static void main(String[] args)
    {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
              
              String pass;
 
              Scanner in = new Scanner(System.in);
              
              System.out.println("Enter passwd");
              pass = in.nextLine();
              System.out.println("You entered "+pass);
              
              try {
                client omeroclient = new client("omero.bioinformatics.ic.ac.uk", 4064);
                //client omeroclient = new client("localhost", 4064);
                //client omeroclient = new client("demo.openmicroscopy.org", 4064);
                

                ServiceFactoryPrx session = omeroclient.createSession("imunro", pass);
                
                long uId = session.getAdminService().getEventContext().userId;  
               
                if (omeroclient != null)  {
                
                  //OMEROImageChooser chooser = new OMEROImageChooser(omeroclient, uId, new Long(4477));
                  String[] strings = {"","filename",".xml"};
                  int type = 1;
                  OMEROImageChooser chooser = new OMEROImageChooser(omeroclient, uId , new  Long(51), strings);
                  
                
                 Dataset returned = chooser.getSelectedDataset();
                // Plate returned = chooser.getSelectedPlate();
                //  OriginalFile returnedFile = chooser.getSelectedFile();
                 
                 if (returned != null)  {
                  System.out.println(returned.getName().getValue());
                  System.out.println(chooser.getFilename() );
                 } 
                 
                 /*Image[] returned = chooser.getSelectedImages();
                  for (int i = 0; i < returned.length; i++) {
                    System.out.println(returned[i].getName().getValue());
                  }  */
                
                  System.out.println("closing down");
                     

                  omeroclient.closeSession();
                }

              } catch (CannotCreateSessionException ex) {
                Logger.getLogger(OMEROImageChooser.class.getName()).log(Level.SEVERE, null, ex);
              } catch (PermissionDeniedException ex) {
                Logger.getLogger(OMEROImageChooser.class.getName()).log(Level.SEVERE, null, ex);
              } catch (ServerError ex) {
                Logger.getLogger(OMEROImageChooser.class.getName()).log(Level.SEVERE, null, ex);
              }
            }

        });
    }
}

