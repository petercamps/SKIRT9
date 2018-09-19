/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "WizardEngine.hpp"
#include "BasicChoiceWizardPane.hpp"
#include "BoolPropertyHandler.hpp"
#include "BoolPropertyWizardPane.hpp"
#include "CreateRootWizardPane.hpp"
#include "DoubleListPropertyHandler.hpp"
#include "DoubleListPropertyWizardPane.hpp"
#include "DoublePropertyHandler.hpp"
#include "DoublePropertyWizardPane.hpp"
#include "EnumPropertyHandler.hpp"
#include "EnumPropertyWizardPane.hpp"
#include "IntPropertyHandler.hpp"
#include "IntPropertyWizardPane.hpp"
#include "Item.hpp"
#include "ItemListPropertyHandler.hpp"
#include "ItemListPropertyWizardPane.hpp"
#include "ItemPropertyHandler.hpp"
#include "ItemPropertyWizardPane.hpp"
#include "ItemUtils.hpp"
#include "MultiPropertyWizardPane.hpp"
#include "OpenWizardPane.hpp"
#include "PropertyHandlerVisitor.hpp"
#include "SaveWizardPane.hpp"
#include "SchemaDef.hpp"
#include "StringPropertyHandler.hpp"
#include "StringPropertyWizardPane.hpp"
#include "StringUtils.hpp"
#include "SubItemPropertyWizardPane.hpp"

////////////////////////////////////////////////////////////////////

WizardEngine::WizardEngine(QObject* parent)
    : QObject(parent)
{
}

////////////////////////////////////////////////////////////////////

WizardEngine::~WizardEngine()
{
}

////////////////////////////////////////////////////////////////////

bool WizardEngine::canAdvance()
{
    switch (_stage)
    {
    case Stage::BasicChoice:
        return _schema != nullptr;
    case Stage::CreateRoot:
        return (_schema && _root) ? _schema->inherits(_root->type(), _schema->schemaType()) : false;
    case Stage::OpenHierarchy:
        return !_filepath.isEmpty();
    case Stage::ConstructHierarchy:
        return _propertyValid;
    case Stage::SaveHierarchy:
        return false;
    }
    return false;   // to satisfy some compilers
}

////////////////////////////////////////////////////////////////////

bool WizardEngine::canRetreat()
{
    return _stage != Stage::BasicChoice;
}

////////////////////////////////////////////////////////////////////

bool WizardEngine::isDirty()
{
    return _dirty;
}

////////////////////////////////////////////////////////////////////

QString WizardEngine::filepath()
{
    return _filepath;
}

////////////////////////////////////////////////////////////////////

const SchemaDef* WizardEngine::schema()
{
    return _schema.get();
}

////////////////////////////////////////////////////////////////////

std::unique_ptr<PropertyHandler> WizardEngine::createPropertyHandler(int propertyIndex)
{
    return _schema->createPropertyHandler(_current, _schema->properties(_current->type())[propertyIndex], &_nameMgr);
}

////////////////////////////////////////////////////////////////////

int WizardEngine::propertyIndexForChild(Item* child)
{
    int index = 0;
    Item* parent = child->parent();
    if (parent) for (auto property : _schema->properties(parent->type()))
    {
        auto handler = _schema->createPropertyHandler(parent, property, &_nameMgr);

        // check the value of item properties
        auto itemhandler = dynamic_cast<ItemPropertyHandler*>(handler.get());
        if (itemhandler && itemhandler->value() == child) return index;

        // check the values of item list properties
        auto itemlisthandler = dynamic_cast<ItemListPropertyHandler*>(handler.get());
        if (itemlisthandler)
        {
            for (auto item : itemlisthandler->value()) if (item == child) return index;
        }

        index++;
    }
    return -1;  // this should never happen
}

////////////////////////////////////////////////////////////////////

bool WizardEngine::isPropertyEligableForMultiPane(int propertyIndex)
{
    auto handler = createPropertyHandler(propertyIndex);
    if (dynamic_cast<StringPropertyHandler*>(handler.get())) return true;
    if (dynamic_cast<BoolPropertyHandler*>(handler.get())) return true;
    if (dynamic_cast<IntPropertyHandler*>(handler.get())) return true;
    if (dynamic_cast<EnumPropertyHandler*>(handler.get())) return true;
    if (dynamic_cast<DoublePropertyHandler*>(handler.get())) return true;
    if (dynamic_cast<DoubleListPropertyHandler*>(handler.get())) return true;
    return false;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // The functions in this class are part of the visitor pattern initiated by the setupProperties() function.
    // They set the value of a non-compound property to its default value, the value of a compound property
    // to "not present", i.e. the null pointer or the empty list, and the value of an item property that
    // offers only one choice to the "forced" value
    class SilentPropertySetter : public PropertyHandlerVisitor
    {
    public:
        SilentPropertySetter() { }

        void visitPropertyHandler(StringPropertyHandler* handler) override
        {
            if (handler->hasDefaultValue())
            {
                handler->setValue(handler->defaultValue());
                handler->setConfigured();
            }
            else handler->setValue(string());
        }

        void visitPropertyHandler(BoolPropertyHandler* handler) override
        {
            if (handler->hasDefaultValue())
            {
                handler->setValue(handler->defaultValue());
                handler->setConfigured();
            }
            else handler->setValue(false);
        }

        void visitPropertyHandler(IntPropertyHandler* handler) override
        {
            if (handler->hasDefaultValue())
            {
                handler->setValue(handler->defaultValue());
                handler->setConfigured();
            }
            else handler->setValue(0);
        }

        void visitPropertyHandler(EnumPropertyHandler* handler) override
        {
            if (handler->hasDefaultValue())
            {
                handler->setValue(handler->defaultValue());
                handler->setConfigured();
            }
            else handler->setValue(handler->values()[0]);
        }

        void visitPropertyHandler(DoublePropertyHandler* handler) override
        {
            if (handler->hasDefaultValue())
            {
                handler->setValue(handler->defaultValue());
                handler->setConfigured();
            }
            else handler->setValue(0.);
        }

        void visitPropertyHandler(DoubleListPropertyHandler* handler) override
        {
            if (handler->hasDefaultValue())
            {
                handler->setValue(handler->defaultValue());
                handler->setConfigured();
            }
            else handler->setValue(vector<double>());
        }

        void visitPropertyHandler(ItemPropertyHandler* handler) override
        {
            if (handler->isRelevant())
            {
                if (handler->hasDefaultValue())
                {
                    handler->setToNewItemOfType(handler->defaultType());
                    handler->setConfigured();
                    return;
                }
                if (handler->isRequired())
                {
                    auto choices = handler->allowedAndDisplayedDescendants();
                    if (choices.size() == 1)
                    {
                        handler->setToNewItemOfType(choices[0]);
                        handler->setConfigured();
                        return;
                    }
                }
            }
            handler->setToNull();
        }

        void visitPropertyHandler(ItemListPropertyHandler* handler) override
        {
            if (handler->isRelevant() && handler->isRequired())
            {
                if (handler->hasDefaultValue())
                {
                    handler->addNewItemOfType(handler->defaultType());
                    handler->setConfigured();
                    return;
                }
                auto choices = handler->allowedAndDisplayedDescendants();
                if (choices.size() == 1)
                {
                    handler->addNewItemOfType(choices[0]);
                    handler->setConfigured();
                    return;
                }
            }
            handler->setToEmpty();
        }
    };
}

////////////////////////////////////////////////////////////////////

bool WizardEngine::isPropertySilent(PropertyHandler* handler)
{
    // an irrelevant property is always silent
    if (!handler->isRelevant()) return true;

    // a property that should not be displayed is silent unless it is required and has no default value
    if (!handler->isDisplayed() && (!handler->isRequired() || handler->hasDefaultValue())) return true;

    // an item property that offers only a single choice is silent
    auto itemhdlr = dynamic_cast<ItemPropertyHandler*>(handler);
    if (itemhdlr)
    {
        auto numChoices = itemhdlr->allowedAndDisplayedDescendants().size();
        if (numChoices == 0) return true;
        if (numChoices == 1 && itemhdlr->isRequired()) return true;
    }

    // the subitem for an item list property that offers only a single choice is silent
    auto itemlisthdlr = dynamic_cast<ItemListPropertyHandler*>(handler);
    if (itemlisthdlr && _subItemIndex < 0)
    {
        if (itemlisthdlr->allowedAndDisplayedDescendants().size() <= 1) return true;
    }

    // if we reach here, the property is not silent
    return false;
}

////////////////////////////////////////////////////////////////////

bool WizardEngine::isCurrentPropertyRangeSilent()
{
    // initialize name manager up to just before the current property range
    createPropertyHandler(_firstPropertyIndex)->rebuildNames();

    // becomes false if at least one of the properties in the range is not silent
    bool result = true;

    // loop over all properties in the range
    for (int propertyIndex=_firstPropertyIndex; propertyIndex<=_lastPropertyIndex; ++propertyIndex)
    {
        auto handler = createPropertyHandler(propertyIndex);

        // if this property is not silent, the complete range is not silent
        if (!isPropertySilent(handler.get())) result = false;

        // if this silent property was not previously configured by the user, set its default value;
        // the corresponding names are automatically inserted
        if (!handler->isConfigured())
        {
            SilentPropertySetter silentSetter;
            handler->acceptVisitor(&silentSetter);
        }
        // otherwise, explicitly insert the names for the property
        else handler->insertNames();
    }
    return result;
}

////////////////////////////////////////////////////////////////////

void WizardEngine::advance(bool recursive)
{
    // remember the current state so we can retreat to it
    if (!recursive) _stateStack.emplace(_stage, _current, _firstPropertyIndex, _lastPropertyIndex, _subItemIndex);

    // advance the state depending on the current stage and details within the stage
    switch (_stage)
    {
    case Stage::BasicChoice:
        {
            _stage = _openExisting ? Stage::OpenHierarchy : Stage::CreateRoot;
            break;
        }
    case Stage::OpenHierarchy:
        {
            _stage = Stage::CreateRoot;
            break;
        }
    case Stage::CreateRoot:
        {
            _stage = Stage::ConstructHierarchy;
            _current = _root.get();
            _firstPropertyIndex = 0; // assumes that the root has at least one property
            break;
        }
    case Stage::ConstructHierarchy:
        {
            // if the (single) property being handled is an item or an item list, we may need to descend the hierarchy
            if (_lastPropertyIndex == _firstPropertyIndex)
            {
                auto handler = createPropertyHandler(_firstPropertyIndex);
                auto itemhdlr = dynamic_cast<ItemPropertyHandler*>(handler.get());
                auto itemlisthdlr = dynamic_cast<ItemListPropertyHandler*>(handler.get());

                // if the property being handled is an item, and the item has properties, then descend the hierarchy
                if (itemhdlr && itemhdlr->value() && _schema->properties(itemhdlr->value()->type()).size()>0)
                {
                    _current = itemhdlr->value();
                    _firstPropertyIndex = 0;
                    break;
                }

                // if the property being handled is an item list, and we're editing one of its subitems,
                // and the subitem has properties, then descend the hierarchy into that subitem
                if (itemlisthdlr && _subItemIndex>=0 &&
                    _schema->properties(itemlisthdlr->value()[_subItemIndex]->type()).size()>0)
                {
                    _current = itemlisthdlr->value()[_subItemIndex];
                    _firstPropertyIndex = 0;
                    break;
                }
            }

            // if we did not descend the hierarchy, attempt to advance to the next property
            _firstPropertyIndex = _lastPropertyIndex+1;

            // if we handled the last property at this level, move up the hierarchy to a level where
            // there are properties to advance to; if we encounter the root item, then move to the SaveHierarchy stage
            while (static_cast<size_t>(_firstPropertyIndex) == _schema->properties(_current->type()).size())
            {
                // indicate that the item we're backing out of is "complete"
                ItemUtils::setItemComplete(_current);

                // special case for root
                if (_current == _root.get())
                {
                    _stage = Stage::SaveHierarchy;
                    break;
                }

                // move up the hierarchy
                _firstPropertyIndex = propertyIndexForChild(_current);
                _current = _current->parent();

                // if we're advancing out of a subitem, stay with the item list property
                if (dynamic_cast<ItemListPropertyHandler*>(createPropertyHandler(_firstPropertyIndex).get()))
                {
                    // and also chop off the retreat states for the sub-item editing sequence
                    while (_stateStack.size() > _stateIndexStack.top()) _stateStack.pop();
                    _stateIndexStack.pop();
                }
                // otherwise go to the next property
                else _firstPropertyIndex++;
            }
            break;
        }
    case Stage::SaveHierarchy:
        {
            break;
        }
    }

    if (_stage == Stage::CreateRoot)
    {
        // skip create root pane if it offers only one choice
        auto choices = _schema->descendants(_schema->schemaType());
        if (choices.size() == 1)
        {
            if (!_root || !_schema->inherits(_root->type(), _schema->schemaType()))
            {
                setRootType(choices[0]);
            }
            advance(true);
        }
    }
    else if (_stage == Stage::ConstructHierarchy)
    {
        // a regular advance can never descend into a subitem, so we always clear the sub-item index
        // (this is meaningless and harmless if the current item is not an item list)
        _subItemIndex = -1;

        // determine the range of properties that can be combined onto a single multi-pane
        _lastPropertyIndex = _firstPropertyIndex;
        if (isPropertyEligableForMultiPane(_firstPropertyIndex))
        {
            while (static_cast<size_t>(_lastPropertyIndex+1) != _schema->properties(_current->type()).size()
                   && isPropertyEligableForMultiPane(_lastPropertyIndex+1)) _lastPropertyIndex++;
        }

        // skip silent properties after setting their default values
        if (isCurrentPropertyRangeSilent())
        {
            advance(true);
        }
    }

    if (!recursive) emitStateChanged();
}


////////////////////////////////////////////////////////////////////

void WizardEngine::advanceToEditSubItem(int subItemIndex)
{
    // remember the index of the previous state,
    // so that we can chop off all subsequent states when advancing out of a sub-item
    // (with the effect that subsequent retreats do not goi back into the sub-item editing sequence)
    _stateIndexStack.emplace(_stateStack.size());

    // remember the current state so we can retreat to it
    _stateStack.emplace(_stage, _current, _firstPropertyIndex, _lastPropertyIndex, _subItemIndex);

    // indicate that we're editing the current sub-item
    _subItemIndex = subItemIndex;

    // skip this wizard pane if there is only one choice for the subitem class
    if (isCurrentPropertyRangeSilent()) advance(true);

    emitStateChanged();
}

////////////////////////////////////////////////////////////////////

void WizardEngine::retreat()
{
    // restore the previous state
    _stateStack.top().getState(_stage, _current, _firstPropertyIndex, _lastPropertyIndex, _subItemIndex);
    _stateStack.pop();

    emitStateChanged();
}

////////////////////////////////////////////////////////////////////

void WizardEngine::emitStateChanged()
{
    emit stateChanged();
    emit canAdvanceChangedTo(canAdvance());
    emit canRetreatChangedTo(canRetreat());
}

////////////////////////////////////////////////////////////////////

void WizardEngine::setBasicChoice(bool openExisting, string libraryPath, string schemaName)
{
    if (_openExisting != openExisting || _schemaName != schemaName)
    {
        // update the choice
        _openExisting = openExisting;
        _schemaName = schemaName;
        _schema = std::make_unique<SchemaDef>(StringUtils::joinPaths(libraryPath, schemaName));

        // clear the current hierarchy and the related state
        _root.reset();
        _filepath.clear();
        _dirty = false;
        emit titleChanged();
        emit dirtyChanged();
        emit canAdvanceChangedTo(canAdvance());
    }
}

////////////////////////////////////////////////////////////////////

void WizardEngine::setRootType(string newRootType)
{
    if (_root && _root->type()==newRootType) return;
    _root = _schema->createItem(newRootType);
    emit canAdvanceChangedTo(canAdvance());
    _dirty = true;
    emit dirtyChanged();
}

////////////////////////////////////////////////////////////////////

void WizardEngine::hierarchyWasLoaded(Item* root, QString filepath)
{
    _root.reset(root);
    hierarchyWasSaved(filepath);
}

////////////////////////////////////////////////////////////////////

void WizardEngine::setPropertyValid(bool valid)
{
    _propertyValid = valid;
    emit canAdvanceChangedTo(canAdvance());
}

////////////////////////////////////////////////////////////////////

void WizardEngine::hierarchyWasChanged()
{
    _dirty = true;
    emit dirtyChanged();
    ItemUtils::setItemIncomplete(_current);
}

////////////////////////////////////////////////////////////////////

void WizardEngine::hierarchyWasSaved(QString filepath)
{
    _filepath = filepath;
    _dirty = false;
    emit titleChanged();
    emit dirtyChanged();
    emit canAdvanceChangedTo(canAdvance());
}

////////////////////////////////////////////////////////////////////

void WizardEngine::restartWizard()
{
    // clear the state
    _openExisting = false;
    _schemaName.clear();
    _schema.reset();
    _root.reset();
    _stage = Stage::BasicChoice;
    _current = nullptr;
    while (!_stateStack.empty()) _stateStack.pop();
    while (!_stateIndexStack.empty()) _stateIndexStack.pop();
    _dirty = false;
    _filepath.clear();
    _nameMgr.clearAll();

    // emit notifications
    emit titleChanged();
    emit dirtyChanged();
    emitStateChanged();
}

////////////////////////////////////////////////////////////////////

QWidget* WizardEngine::createPane()
{
    switch (_stage)
    {
    case Stage::BasicChoice:
        {
            return new BasicChoiceWizardPane(_openExisting, _schemaName, _dirty, this);
        }
    case Stage::CreateRoot:
        {
            string currentType = _root ? _root->type() : "";
            return new CreateRootWizardPane(_schema.get(), currentType, this);
        }
    case Stage::OpenHierarchy:
        {
            return new OpenWizardPane(_schema.get(), _filepath, _dirty, this);
        }
    case Stage::ConstructHierarchy:
        {
            // single pane
            if (_lastPropertyIndex == _firstPropertyIndex)
            {
                auto handler = createPropertyHandler(_firstPropertyIndex);

                if (dynamic_cast<StringPropertyHandler*>(handler.get()))
                    return new StringPropertyWizardPane(std::move(handler), this);
                if (dynamic_cast<BoolPropertyHandler*>(handler.get()))
                    return new BoolPropertyWizardPane(std::move(handler), this);
                if (dynamic_cast<IntPropertyHandler*>(handler.get()))
                    return new IntPropertyWizardPane(std::move(handler), this);
                if (dynamic_cast<EnumPropertyHandler*>(handler.get()))
                    return new EnumPropertyWizardPane(std::move(handler), this);
                if (dynamic_cast<DoublePropertyHandler*>(handler.get()))
                    return new DoublePropertyWizardPane(std::move(handler), this);
                if (dynamic_cast<DoubleListPropertyHandler*>(handler.get()))
                    return new DoubleListPropertyWizardPane(std::move(handler), this);
                if (dynamic_cast<ItemPropertyHandler*>(handler.get()))
                    return new ItemPropertyWizardPane(std::move(handler), this);
                if (dynamic_cast<ItemListPropertyHandler*>(handler.get()))
                {
                    if (_subItemIndex<0) return new ItemListPropertyWizardPane(std::move(handler), this);
                    else                 return new SubItemPropertyWizardPane(std::move(handler), this);
                }
            }

            // multi-pane
            else
            {
                auto multipane = new MultiPropertyWizardPane(this);
                for (int propertyIndex=_firstPropertyIndex; propertyIndex<=_lastPropertyIndex; ++propertyIndex)
                {
                    auto handler = createPropertyHandler(propertyIndex);

                    if (dynamic_cast<StringPropertyHandler*>(handler.get()))
                        multipane->addPane(new StringPropertyWizardPane(std::move(handler), multipane));
                    else if (dynamic_cast<BoolPropertyHandler*>(handler.get()))
                        multipane->addPane(new BoolPropertyWizardPane(std::move(handler), multipane));
                    else if (dynamic_cast<IntPropertyHandler*>(handler.get()))
                        multipane->addPane(new IntPropertyWizardPane(std::move(handler), multipane));
                    else if (dynamic_cast<EnumPropertyHandler*>(handler.get()))
                        multipane->addPane(new EnumPropertyWizardPane(std::move(handler), multipane));
                    else if (dynamic_cast<DoublePropertyHandler*>(handler.get()))
                        multipane->addPane(new DoublePropertyWizardPane(std::move(handler), multipane));
                    else if (dynamic_cast<DoubleListPropertyHandler*>(handler.get()))
                        multipane->addPane(new DoubleListPropertyWizardPane(std::move(handler), multipane));
                }
                return multipane;
            }
            break;  // to satisfy gcc compiler
        }
    case Stage::SaveHierarchy:
        {
            return new SaveWizardPane(_schema.get(), _root.get(), _filepath, _dirty, this);
        }
    }
    return nullptr;
}

////////////////////////////////////////////////////////////////////

QString WizardEngine::hierarchyPath()
{
    string result;

    if (_stage == Stage::ConstructHierarchy)
    {
        // on the lowest level, show item type and property name
        result = _current->type() + " : " + _schema->properties(_current->type())[_firstPropertyIndex];

        // for higher levels, show only item type
        Item* current = _current->parent();
        while (current)
        {
            result = current->type() + u8" \u2192 " + result;
            current = current->parent();
        }
    }
    return QString::fromStdString(result);
}

////////////////////////////////////////////////////////////////////
