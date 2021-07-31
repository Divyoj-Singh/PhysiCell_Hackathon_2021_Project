/**
 * @file   Output.h
 * @brief  Implementation of the Output class
 * @author Generated by autocreate code
 *
 * <!--------------------------------------------------------------------------
 * This file is part of libSBML.  Please visit http://sbml.org for more
 * information about SBML, and the latest version of libSBML.
 *
 * Copyright (C) 2013-2016 jointly by the following organizations:
 *     1. California Institute of Technology, Pasadena, CA, USA
 *     2. EMBL European Bioinformatics Institute (EMBL-EBI), Hinxton, UK
 *     3. University of Heidelberg, Heidelberg, Germany
 * 
 * Copyright (C) 2009-2013 jointly by the following organizations:
 *     1. California Institute of Technology, Pasadena, CA, USA
 *     2. EMBL European Bioinformatics Institute (EMBL-EBI), Hinxton, UK
 *
 * Copyright (C) 2006-2008 by the California Institute of Technology,
 *     Pasadena, CA, USA 
 *
 * Copyright (C) 2002-2005 jointly by the following organizations:
 *     1. California Institute of Technology, Pasadena, CA, USA
 *     2. Japan Science and Technology Agency, Japan
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution
 * and also available online as http://sbml.org/software/libsbml/license.html
 * ------------------------------------------------------------------------ -->
 *
 * @class Output
 * @sbmlbrief{qual} The output from a QualitativeSpecies.
 * 
 * Each Output refers to a QualitativeSpecies that participates in (is
 * affected by) the corresponding Transition. In Petri net models these are
 * the output places of the transition.
 *
 * In a logical model, a QualitativeSpecies should be referenced in at most
 * one ListOfOutputs, (that of the Transition defining the evolution of this
 * species). When a Transition has several outputs, it is because the
 * referenced species share the same regulators and the same logical rules.
 */


#ifndef Output_H__
#define Output_H__


#include <sbml/common/extern.h>
#include <sbml/common/sbmlfwd.h>
#include <sbml/packages/qual/common/qualfwd.h>

LIBSBML_CPP_NAMESPACE_BEGIN

/** 
 * @enum OutputTransitionEffect_t
 * @brief Enumeration of possible values for the 'transitionEffect' attribute of an Output.
 */
typedef enum
{
    OUTPUT_TRANSITION_EFFECT_PRODUCTION /*!< 'production':  The level of the QualitativeSpecies is increased by the resultLevel of the applicable FunctionTerm possibly modified by the outputLevel of the Output. */
  , OUTPUT_TRANSITION_EFFECT_ASSIGNMENT_LEVEL /*!< 'assignmentLevel':  The level of the qualitativeSpecies is set to the resultLevel of the selected term. */
  , OUTPUT_TRANSITION_EFFECT_UNKNOWN /*!< Unknown or illegal value:  anything other than 'production' or 'assignmentLevel'. */
} OutputTransitionEffect_t;

LIBSBML_CPP_NAMESPACE_END


#ifdef __cplusplus


#include <string>


#include <sbml/SBase.h>
#include <sbml/ListOf.h>
#include <sbml/packages/qual/extension/QualExtension.h>


LIBSBML_CPP_NAMESPACE_BEGIN


class LIBSBML_EXTERN Output : public SBase
{

protected:

  /** @cond doxygenLibsbmlInternal */
//  std::string   mId;
  std::string   mQualitativeSpecies;
  OutputTransitionEffect_t   mTransitionEffect;
//  std::string   mName;
  int           mOutputLevel;
  bool          mIsSetOutputLevel;
  /** @endcond */


public:

  /**
   * Creates a new Output with the given level, version, and package version.
   *
   * @param level an unsigned int, the SBML Level to assign to this Output.
   *
   * @param version an unsigned int, the SBML Version to assign to this Output.
   *
   * @param pkgVersion an unsigned int, the SBML Qual Version to assign to this Output.
   *
   * @copydetails doc_note_setting_lv_pkg
   */
  Output(unsigned int level      = QualExtension::getDefaultLevel(),
         unsigned int version    = QualExtension::getDefaultVersion(),
         unsigned int pkgVersion = QualExtension::getDefaultPackageVersion());


  /**
   * Creates a new Output with the given QualPkgNamespaces object.
   *
   * @copydetails doc_what_are_sbml_package_namespaces
   *
   * @param qualns the QualPkgNamespaces object.
   *
   * @copydetails doc_note_setting_lv_pkg
   */
  Output(QualPkgNamespaces* qualns);


  /**
   * Copy constructor for Output.
   *
   * @param orig the Output instance to copy.
   */
  Output(const Output& orig);


  /**
   * Assignment operator for Output.
   *
   * @param rhs the object whose values are used as the basis
   * of the assignment.
   */
  Output& operator=(const Output& rhs);


  /**
   * Creates and returns a deep copy of this Output object.
   *
   * @return a (deep) copy of this Output object.
   */
  virtual Output* clone () const;


  /**
   * Destructor for Output.
   */
  virtual ~Output();


  /**
   * Returns the value of the "id" attribute of this Output.
   *
   * @note Because of the inconsistent behavior of this function with 
   * respect to assignments and rules, it is now recommended to
   * use the getIdAttribute() function instead.
   *
   * @copydetails doc_id_attribute
   *
   * @return the id of this Output.
   *
   * @see getIdAttribute()
   * @see setIdAttribute(const std::string& sid)
   * @see isSetIdAttribute()
   * @see unsetIdAttribute()
   */
  virtual const std::string& getId() const;


  /**
   * Returns the value of the "qualitativeSpecies" attribute of this Output.
   *
   * @return the value of the "qualitativeSpecies" attribute of this Output as a string.
   */
  virtual const std::string& getQualitativeSpecies() const;


  /**
   * Returns the value of the "transitionEffect" attribute of this Output.
   *
   * @return the value of the "transitionEffect" attribute of this Output as a string.
   */
  OutputTransitionEffect_t getTransitionEffect() const;


  /**
   * Returns the value of the "name" attribute of this Output object.
   *
   * @copydetails doc_get_name
   */
  virtual const std::string& getName() const;


  /**
   * Returns the value of the "outputLevel" attribute of this Output.
   *
   * @return the value of the "outputLevel" attribute of this Output as a integer.
   */
  virtual int getOutputLevel() const;


  /**
   * Predicate returning @c true or @c false depending on whether this
   * Output's "id" attribute has been set.
   *
   * @copydetails doc_isset_id
   */
  virtual bool isSetId() const;


  /**
   * Predicate returning @c true or @c false depending on whether this
   * Output's "qualitativeSpecies" attribute has been set.
   *
   * @return @c true if this Output's "qualitativeSpecies" attribute has been set,
   * otherwise @c false is returned.
   */
  virtual bool isSetQualitativeSpecies() const;


  /**
   * Predicate returning @c true or @c false depending on whether this
   * Output's "transitionEffect" attribute has been set.
   *
   * @return @c true if this Output's "transitionEffect" attribute has been set,
   * otherwise @c false is returned.
   */
  bool isSetTransitionEffect() const;


  /**
   * Predicate returning @c true or @c false depending on whether this
   * Output's "name" attribute has been set.
   *
   * @copydetails doc_isset_name
   */
  virtual bool isSetName() const;


  /**
   * Predicate returning @c true or @c false depending on whether this
   * Output's "outputLevel" attribute has been set.
   *
   * @return @c true if this Output's "outputLevel" attribute has been set,
   * otherwise @c false is returned.
   */
  virtual bool isSetOutputLevel() const;


  /**
   * Sets the value of the "id" attribute of this Output.
   *
   * @copydetails doc_set_id
   */
  virtual int setId(const std::string& sid);


  /**
   * Sets the value of the "qualitativeSpecies" attribute of this Output.
   *
   * @param qualitativeSpecies the value of the "qualitativeSpecies" attribute to be set.
   *
   * @copydetails doc_returns_success_code
   * @li @sbmlconstant{LIBSBML_OPERATION_SUCCESS, OperationReturnValues_t}
   * @li @sbmlconstant{LIBSBML_INVALID_ATTRIBUTE_VALUE, OperationReturnValues_t}
   */
  virtual int setQualitativeSpecies(const std::string& qualitativeSpecies);


  /**
   * Sets the value of the "transitionEffect" attribute of this Output.
   *
   * @param transitionEffect the value of the "transitionEffect" attribute to be set.
   *
   * @copydetails doc_returns_success_code
   * @li @sbmlconstant{LIBSBML_OPERATION_SUCCESS, OperationReturnValues_t}
   * @li @sbmlconstant{LIBSBML_INVALID_ATTRIBUTE_VALUE, OperationReturnValues_t}
   */
  int setTransitionEffect(const OutputTransitionEffect_t transitionEffect);


  /**
   * Sets the value of the "name" attribute of this Output.
   *
   * @copydetails doc_set_name
   */
  virtual int setName(const std::string& name);


  /**
   * Sets the value of the "outputLevel" attribute of this Output.
   *
   * @param outputLevel int value of the "outputLevel" attribute to be set.
   *
   * @copydetails doc_returns_success_code
   * @li @sbmlconstant{LIBSBML_OPERATION_SUCCESS, OperationReturnValues_t}
   * @li @sbmlconstant{LIBSBML_INVALID_ATTRIBUTE_VALUE, OperationReturnValues_t}
   */
  virtual int setOutputLevel(int outputLevel);


  /**
   * Unsets the value of the "id" attribute of this Output.
   *
   * @copydetails doc_unset_id
   */
  virtual int unsetId();


  /**
   * Unsets the value of the "qualitativeSpecies" attribute of this Output.
   *
   * @copydetails doc_returns_success_code
   * @li @sbmlconstant{LIBSBML_OPERATION_SUCCESS, OperationReturnValues_t}
   * @li @sbmlconstant{LIBSBML_OPERATION_FAILED, OperationReturnValues_t}
   */
  virtual int unsetQualitativeSpecies();


  /**
   * Unsets the value of the "transitionEffect" attribute of this Output.
   *
   * @copydetails doc_returns_success_code
   * @li @sbmlconstant{LIBSBML_OPERATION_SUCCESS, OperationReturnValues_t}
   * @li @sbmlconstant{LIBSBML_OPERATION_FAILED, OperationReturnValues_t}
   */
  int unsetTransitionEffect();


  /**
   * Unsets the value of the "name" attribute of this Output.
   *
   * @copydetails doc_unset_name
   */
  virtual int unsetName();


  /**
   * Unsets the value of the "outputLevel" attribute of this Output.
   *
   * @copydetails doc_returns_success_code
   * @li @sbmlconstant{LIBSBML_OPERATION_SUCCESS, OperationReturnValues_t}
   * @li @sbmlconstant{LIBSBML_OPERATION_FAILED, OperationReturnValues_t}
   */
  virtual int unsetOutputLevel();


  /**
   * @copydoc doc_renamesidref_common
   */
   virtual void renameSIdRefs(const std::string& oldid, const std::string& newid);


  /**
   * Returns the XML name of this object.
   *
   * @return the name of this element.
   */
  virtual const std::string& getElementName () const;


  /**
   * Returns the libSBML type code of this object instance.
   *
   * @copydetails doc_what_are_typecodes
   *
   * @return the SBML type code for this object:
   * @sbmlconstant{SBML_QUAL_OUTPUT, SBMLQualTypeCode_t}
   *
   * @copydetails doc_warning_typecodes_not_unique
   *
   * @see getElementName()
   * @see getPackageName()
   */
  virtual int getTypeCode () const;


  /**
   * Predicate returning @c true if all the required attributes
   * for this Output object have been set.
   *
   * @note The required attributes for a Output object are:
   * @li "qualitativeSpecies"
   * @li "transitionEffect"
   *
   * @return a boolean value indicating whether all the required
   * attributes for this object have been defined.
   */
  virtual bool hasRequiredAttributes() const;


  /** @cond doxygenLibsbmlInternal */
  /**
   * Subclasses should override this method to write out their contained
   * SBML objects as XML elements.  Be sure to call your parent's
   * implementation of this method as well.
   */
  virtual void writeElements (XMLOutputStream& stream) const;


  /** @endcond */


  /** @cond doxygenLibsbmlInternal */
  /**
   * Accepts the given SBMLVisitor.
   */
  virtual bool accept (SBMLVisitor& v) const;
  /** @endcond */


  /** @cond doxygenLibsbmlInternal */
  /**
   * Sets the parent SBMLDocument.
   */
  virtual void setSBMLDocument (SBMLDocument* d);


  /** @endcond */


  /** @cond doxygenLibsbmlInternal */
  /**
   * Enables/Disables the given package with this element.
   */
  virtual void enablePackageInternal(const std::string& pkgURI,
               const std::string& pkgPrefix, bool flag);


  /** @endcond */


protected:

  /** @cond doxygenLibsbmlInternal */
  /**
   * Get the list of expected attributes for this element.
   */
  virtual void addExpectedAttributes(ExpectedAttributes& attributes);


  /** @endcond */


  /** @cond doxygenLibsbmlInternal */
  /**
   * Read values from the given XMLAttributes set into their specific fields.
   */
  virtual void readAttributes (const XMLAttributes& attributes,
                               const ExpectedAttributes& expectedAttributes);


  /** @endcond */


  /** @cond doxygenLibsbmlInternal */
  /**
   * Write values of XMLAttributes to the output stream.
   */
  virtual void writeAttributes (XMLOutputStream& stream) const;


  /** @endcond */



};

/**
 * @class ListOfOutputs
 * @sbmlbrief{qual} A list of Output objects.
 * 
 * The ListOfOutputs is a container for the Output elements of a Transition.
 * 
 * @copydetails doc_what_is_listof
 *
 * @see Output
 */
class LIBSBML_EXTERN ListOfOutputs : public ListOf
{

public:

  /**
   * Creates a new ListOfOutputs with the given level, version, and package version.
   *
   * @param level an unsigned int, the SBML Level to assign to this ListOfOutputs.
   *
   * @param version an unsigned int, the SBML Version to assign to this ListOfOutputs.
   *
   * @param pkgVersion an unsigned int, the SBML Qual Version to assign to this ListOfOutputs.
   *
   * @copydetails doc_note_setting_lv_pkg
   */
  ListOfOutputs(unsigned int level      = QualExtension::getDefaultLevel(),
                unsigned int version    = QualExtension::getDefaultVersion(),
                unsigned int pkgVersion = QualExtension::getDefaultPackageVersion());


  /**
   * Creates a new ListOfOutputs with the given QualPkgNamespaces object.
   *
   * @copydetails doc_what_are_sbml_package_namespaces
   *
   * @param qualns the QualPkgNamespaces object.
   *
   * @copydetails doc_note_setting_lv_pkg
   */
  ListOfOutputs(QualPkgNamespaces* qualns);


  /**
   * Creates and returns a deep copy of this ListOfOutputs object.
   *
   * @return a (deep) copy of this ListOfOutputs object.
   */
  virtual ListOfOutputs* clone () const;


  /**
   * Get a Output from the ListOfOutputs.
   *
   * @param n the index number of the Output to get.
   *
   * @return the nth Output in this ListOfOutputs.
   *
   * @see size()
   */
  virtual Output* get(unsigned int n);


  /**
   * Get a Output from the ListOfOutputs.
   *
   * @param n the index number of the Output to get.
   *
   * @return the nth Output in this ListOfOutputs.
   *
   * @see size()
   */
  virtual const Output* get(unsigned int n) const;


  /**
   * Get a Output from the ListOfOutputs
   * based on its identifier.
   *
   * @param sid a string representing the identifier
   * of the Output to get.
   *
   * @return Output in this ListOfOutputs
   * with the given id or NULL if no such
   * Output exists.
   *
   * @see get(unsigned int n)   *
   * @see size()
   */
  virtual Output* get(const std::string& sid);


  /**
   * Get a Output from the ListOfOutputs
   * based on its identifier.
   *
   * @param sid a string representing the identifier
   * of the Output to get.
   *
   * @return Output in this ListOfOutputs
   * with the given id or NULL if no such
   * Output exists.
   *
   * @see get(unsigned int n)   *
   * @see size()
   */
  virtual const Output* get(const std::string& sid) const;


  /**
   * Get a Output from the ListOfOutputs
   * based on the qualitativeSpecies to which it refers.
   *
   * @param sid a string representing the qualitativeSpecies attribute
   * of the Output to get.
   *
   * @return the first Output in this ListOfOutputs
   * with the given qualitativeSpecies or NULL if no such
   * Output exists.
   *
   * @see get(unsigned int n)   *
   * @see size()
   */
  const Output* getBySpecies(const std::string& sid) const;


  /**
   * Get a Output from the ListOfOutputs
   * based on the qualitativeSpecies to which it refers.
   *
   * @param sid a string representing the qualitativeSpecies attribute
   * of the Output to get.
   *
   * @return the first Output in this ListOfOutputs
   * with the given qualitativeSpecies or NULL if no such
   * Output exists.
   *
   * @see get(unsigned int n)   *
   * @see size()
   */
  Output* getBySpecies(const std::string& sid);


  /**
   * Removes the nth Output from this ListOfOutputs
   * and returns a pointer to it.
   *
   * The caller owns the returned item and is responsible for deleting it.
   *
   * @param n the index of the Output to remove.
   *
   * @see size()
   */
  virtual Output* remove(unsigned int n);


  /**
   * Removes the Output from this ListOfOutputs with the given identifier
   * and returns a pointer to it.
   *
   * The caller owns the returned item and is responsible for deleting it.
   * If none of the items in this list have the identifier @p sid, then
   * @c NULL is returned.
   *
   * @param sid the identifier of the Output to remove.
   *
   * @return the Output removed. As mentioned above, the caller owns the
   * returned item.
   */
  virtual Output* remove(const std::string& sid);


  /**
   * Returns the XML name of this object.
   *
   * @return the name of this element.
   */
  virtual const std::string& getElementName () const;


  /**
   * Returns the libSBML type code for the SBML objects
   * contained in this ListOf object.
   * 
   * @copydetails doc_what_are_typecodes
   *
   * @return the SBML type code for objects contained in this list:
   * @sbmlconstant{SBML_QUAL_OUTPUT, SBMLTypeCode_t} (default).
   *
   * @see getElementName()
   * @see getPackageName()
   */
  virtual int getItemTypeCode () const;


protected:

  /** @cond doxygenLibsbmlInternal */
  /**
   * Creates a new Output in this ListOfOutputs
   */
  virtual SBase* createObject(XMLInputStream& stream);


  /** @endcond */


  /** @cond doxygenLibsbmlInternal */
  /**
   * Write the namespace for the Qual package.
   */
  virtual void writeXMLNS(XMLOutputStream& stream) const;


  /** @endcond */



};



LIBSBML_CPP_NAMESPACE_END

#endif  /*  __cplusplus  */

#ifndef SWIG

LIBSBML_CPP_NAMESPACE_BEGIN
BEGIN_C_DECLS

/**
 * Creates a new Output_t structure using the given SBML @p level
 * and @p version values.
 *
 * @param level an unsigned int, the SBML Level to assign to this
 * Output_t.
 * @param version an unsigned int, the SBML Version to assign to this
 * Output_t.
 * @param pkgVersion an unsigned int, the SBML 'Qual' package Version to assign to this
 * Output_t.
 *
 * @return a pointer to the newly created Output_t structure.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
Output_t *
Output_create(unsigned int level, unsigned int version,
              unsigned int pkgVersion);


/**
 * Frees the given Output_t structure.
 *
 * @param o the Output_t structure to free.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
void
Output_free(Output_t * o);


/**
 * Returns a copy of the given Output_t structure.
 *
 * @param o the Output_t structure to copy.
 * 
 * @return a (deep) copy of the Output_t.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
Output_t *
Output_clone(Output_t * o);


/**
 * Takes an Output_t structure and returns its identifier.
 *
 * @param o the Output_t structure whose identifier is sought.
 * 
 * @return the identifier of this Output_t, as a pointer to a string.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
char *
Output_getId(Output_t * o);


/**
 * Takes a Output_t structure and returns its qualitativeSpecies.
 *
 * @param o the Output_t whose qualitativeSpecies is sought.
 *
 * @return the qualitativeSpecies of the given Output_t, as a pointer to a string.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
char *
Output_getQualitativeSpecies(Output_t * o);


/**
 * Takes a Output_t structure and returns its transitionEffect as a ENUMTYPE 
 * (an enumerated value of legal possible values for the attribute).
 *
 * @param o the Output_t whose transitionEffect is sought.
 *
 * @return the transitionEffect of the given Output_t, as a ENUMTYPE
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
OutputTransitionEffect_t
Output_getTransitionEffect(Output_t * o);


/**
 * Takes a Output_t structure and returns its name.
 *
 * @param o the Output_t whose name is sought.
 *
 * @return the name of this Output_t, as a pointer to a string.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
char *
Output_getName(Output_t * o);


/**
 * Takes a Output_t structure and returns its outputLevel.
 *
 * @param o the Output_t whose outputLevel is sought.
 *
 * @return the outputLevel attribute of the given Output_t, as an @c int.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_getOutputLevel(Output_t * o);


/**
 * Predicate returning @c true or @c false depending on whether the given
 * Output_t structure's identifier is set.
 *
 * @param o the Output_t structure to query.
 * 
 * @return @c non-zero (true) if the "id" attribute of the given
 * Output_t structure is set, zero (false) otherwise.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_isSetId(Output_t * o);


/**
 * Predicate returning @c true or @c false depending on whether the given
 * Output_t structure's qualitativeSpecies is set.
 *
 * @param o the Output_t structure to query.
 * 
 * @return @c non-zero (true) if the "qualitativeSpecies" attribute of the given
 * Output_t structure is set, zero (false) otherwise.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_isSetQualitativeSpecies(Output_t * o);


/**
 * Predicate returning @c true or @c false depending on whether the given
 * Output_t structure's transitionEffect is set.
 *
 * @param o the Output_t structure to query.
 * 
 * @return @c non-zero (true) if the "transitionEffect" attribute of the given
 * Output_t structure is set, zero (false) otherwise.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_isSetTransitionEffect(Output_t * o);


/**
 * Predicate returning @c true or @c false depending on whether the given
 * Output_t structure's name is set.
 *
 * @param o the Output_t structure to query.
 * 
 * @return @c non-zero (true) if the "name" attribute of the given
 * Output_t structure is set, zero (false) otherwise.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_isSetName(Output_t * o);


/**
 * Predicate returning @c true or @c false depending on whether the given
 * Output_t structure's outputLevel is set.
 *
 * @param o the Output_t structure to query.
 * 
 * @return @c non-zero (true) if the "outputLevel" attribute of the given
 * Output_t structure is set, zero (false) otherwise.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_isSetOutputLevel(Output_t * o);


/**
 * Assigns the identifier of an Output_t structure.
 *
 * This makes a copy of the string passed in the param @p sid.
 *
 * @param o the Output_t structure to set.
 * @param sid the string to use as the identifier.
 *
 * @copydetails doc_returns_success_code
 * @li @sbmlconstant{LIBSBML_OPERATION_SUCCESS, OperationReturnValues_t}
 * @li @sbmlconstant{LIBSBML_INVALID_ATTRIBUTE_VALUE, OperationReturnValues_t}
 *
 * @note Using this function with an id of NULL is equivalent to
 * unsetting the "id" attribute.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_setId(Output_t * o, const char * sid);


/**
 * Sets the qualitativeSpecies of the given Output_t to a copy of @p qualitativeSpecies.
 *
 * @param o the Output_t structure to set.
 * @param qualitativeSpecies the qualitativeSpecies to assign to the given Output_t's "qualitativeSpecies" attribute.
 *
 * @copydetails doc_returns_success_code
 * @li @sbmlconstant{LIBSBML_OPERATION_SUCCESS, OperationReturnValues_t}
 * @li @sbmlconstant{LIBSBML_INVALID_ATTRIBUTE_VALUE, OperationReturnValues_t}
 *
 * @note Using this function with the name set to NULL is equivalent to
 * unsetting the "qualitativeSpecies" attribute.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_setQualitativeSpecies(Output_t * o, const char * qualitativeSpecies);


/**
 * Sets the transitionEffect of the given Output_t to a copy of @p transitionEffect.
 *
 * @param o the Output_t structure to set.
 * @param transitionEffect the transitionEffect to assign to the given Output_t's "transitionEffect" attribute.
 *
 * @copydetails doc_returns_success_code
 * @li @sbmlconstant{LIBSBML_OPERATION_SUCCESS, OperationReturnValues_t}
 * @li @sbmlconstant{LIBSBML_INVALID_ATTRIBUTE_VALUE, OperationReturnValues_t}
 *
 * @note Using this function with the name set to NULL is equivalent to
 * unsetting the "transitionEffect" attribute.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_setTransitionEffect(Output_t * o, const char * transitionEffect);


/**
 * Sets the name of the given Output_t to a copy of @p name.
 *
 * @param o the Output_t structure to set.
 * @param name the name to assign to the given Output_t's "name" attribute.
 *
 * @copydetails doc_returns_success_code
 * @li @sbmlconstant{LIBSBML_OPERATION_SUCCESS, OperationReturnValues_t}
 * @li @sbmlconstant{LIBSBML_INVALID_ATTRIBUTE_VALUE, OperationReturnValues_t}
 *
 * @note Using this function with the name set to NULL is equivalent to
 * unsetting the "name" attribute.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_setName(Output_t * o, const char * name);


/**
 * Sets the "outputLevel" attribute of the given Output_t
 * structure.
 *
 * @param o the Output_t structure.
 * 
 * @param outputLevel the value of outputLevel to assign to the "outputLevel" attribute.
 *
 * @copydetails doc_returns_success_code
 * @li @sbmlconstant{LIBSBML_OPERATION_SUCCESS, OperationReturnValues_t}
 * @li @sbmlconstant{LIBSBML_UNEXPECTED_ATTRIBUTE, OperationReturnValues_t}
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_setOutputLevel(Output_t * o, int outputLevel);


/**
 * Unsets the "id" attribute of the given Output_t structure.
 *
 * @param o the Output_t structure to unset.
 *
 * @copydetails doc_returns_success_code
 * @li @sbmlconstant{LIBSBML_OPERATION_SUCCESS, OperationReturnValues_t}
 * @li @sbmlconstant{LIBSBML_OPERATION_FAILED, OperationReturnValues_t}
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_unsetId(Output_t * o);


/**
 * Unsets the "qualitativeSpecies" attribute of the given Output_t structure.
 *
 * @param o the Output_t structure to unset.
 *
 * @copydetails doc_returns_success_code
 * @li @sbmlconstant{LIBSBML_OPERATION_SUCCESS, OperationReturnValues_t}
 * @li @sbmlconstant{LIBSBML_OPERATION_FAILED, OperationReturnValues_t}
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_unsetQualitativeSpecies(Output_t * o);


/**
 * Unsets the "transitionEffect" attribute of the given Output_t structure.
 *
 * @param o the Output_t structure to unset.
 *
 * @copydetails doc_returns_success_code
 * @li @sbmlconstant{LIBSBML_OPERATION_SUCCESS, OperationReturnValues_t}
 * @li @sbmlconstant{LIBSBML_OPERATION_FAILED, OperationReturnValues_t}
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_unsetTransitionEffect(Output_t * o);


/**
 * Unsets the "name" attribute of the given Output_t structure.
 *
 * @param o the Output_t structure to unset.
 *
 * @copydetails doc_returns_success_code
 * @li @sbmlconstant{LIBSBML_OPERATION_SUCCESS, OperationReturnValues_t}
 * @li @sbmlconstant{LIBSBML_OPERATION_FAILED, OperationReturnValues_t}
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_unsetName(Output_t * o);


/**
 * Unsets the "outputLevel" attribute of the given Output_t structure.
 *
 * @param o the Output_t structure to unset.
 *
 * @copydetails doc_returns_success_code
 * @li @sbmlconstant{LIBSBML_OPERATION_SUCCESS, OperationReturnValues_t}
 * @li @sbmlconstant{LIBSBML_OPERATION_FAILED, OperationReturnValues_t}
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_unsetOutputLevel(Output_t * o);


/**
  * Predicate returning @c true or @c false depending on whether
  * all the required attributes for the given Output_t structure
  * have been set.
  *
  * @note The required attributes for a Output_t structure are:
  * @li useValuesfromTriggerTime ( L3 onwards )
  *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int
Output_hasRequiredAttributes(Output_t * o);


/**
 * Return the Output_t indicated by the given @p sid.
 *
 * @param lo the ListOf_t structure to use.
 *
 * @param sid a string, the identifier of the
 * Output_t is being sought.
 *
 * @return the Output_t for the given variable, or @c NULL if no such
 * Output_t exits.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
Output_t *
ListOfOutputs_getById(ListOf_t * lo, const char * sid);


/**
 * Removes the Output_t structure with the given @p sid
 * from the given ListOf_t structure and returns a pointer to it.
 *
 * The caller owns the returned structure and is responsible for deleting it.
 *
 * @param lo the ListOf_t structure.
 * @param sid the string of the "id" attribute of the Output_t sought.
 *
 * @return the Output_t structure removed.  As mentioned above, the 
 * caller owns the returned structure. @c NULL is returned if no Output_t
 * structure with the "id" attribute exists in the given ListOf_t structure.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
Output_t *
ListOfOutputs_removeById(ListOf_t * lo, const char * sid);

/**
 * Returns the string version of the provided OutputTransitionEffect_t enumeration.
 *
 * @param effect the OutputTransitionEffect_t enumeration to convert.
 *
 * @return A string corresponding to the given effect:  "production", 
 * "assignmentLevel", or @c NULL if the value is OUTPUT_TRANSITION_EFFECT_UNKNOWN 
 * or another invalid enumeration value.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
const char* 
OutputTransitionEffect_toString(OutputTransitionEffect_t effect);


/**
 * Returns the OutputTransitionEffect_t enumeration corresponding to 
 * the given string, or OUTPUT_TRANSITION_EFFECT_UNKNOWN if there is 
 * no such match.  The matching is case-sensitive:  "production" will 
 * return OUTPUT_TRANSITION_EFFECT_PRODUCTION, but "Production" will return 
 * OUTPUT_TRANSITION_EFFECT_UNKNOWN.
 *
 * @param s the string to convert to an OutputTransitionEffect_t.
 *
 * @return The corresponding OutputTransitionEffect_t, or 
 * OUTPUT_TRANSITION_EFFECT_UNKNOWN if no match found.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
OutputTransitionEffect_t 
OutputTransitionEffect_fromString(const char* s);


/**
 * Predicate returning @c true (non-zero) or @c false (zero) depending on whether the given
 * OutputTransitionEffect_t is valid.
 *
 * @param effect the OutputTransitionEffect_t enumeration to query.
 * 
 * @return @c non-zero (true) if the OutputTransitionEffect_t is
 * OUTPUT_TRANSITION_EFFECT_PRODUCTION or OUTPUT_TRANSITION_EFFECT_ASSIGNMENT_LEVEL,
 * zero (false) otherwise (including OUTPUT_TRANSITION_EFFECT_UNKNOWN).
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int 
OutputTransitionEffect_isValidOutputTransitionEffect(OutputTransitionEffect_t effect);


/**
 * Predicate returning @c true (non-zero) or @c false (zero) depending 
 * on whether the given string is a valid OutputTransitionEffect_t.  
 * The matching is case-sensitive:  "production" will return @c true, but 
 * "Production" will return @c false.
 *
 * @param s the string to query.
 * 
 * @return @c non-zero (true) if the string is
 * "production" or "assignmentLevel"; zero (false) otherwise.
 *
 * @memberof Output_t
 */
LIBSBML_EXTERN
int 
OutputTransitionEffect_isValidOutputTransitionEffectString(const char* s);


END_C_DECLS
LIBSBML_CPP_NAMESPACE_END

#endif  /*  !SWIG  */

#endif /*  Output_H__  */
