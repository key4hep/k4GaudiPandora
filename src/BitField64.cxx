#include "BitField64.hxx"

BitFieldValue::BitFieldValue(int64_t& bitfield, const std::string& theName, unsigned theOffset, int signedWidth)
        : _b(bitfield), _mask(0), _name(theName), _offset(theOffset), _width(abs(signedWidth)),
          _minVal(0), _maxVal(0), _isSigned(signedWidth < 0) {

	// Ensure the information fits in the object
        if (_offset > 63 || _offset + _width > 64) {
            throw std::runtime_error("BitFieldValue: offset out of range");
        }
	
	// Setup initial values
        _mask = ((0x0001LL << _width) - 1) << _offset;
	
        if (_isSigned) {
            _minVal = (1LL << (_width - 1)) - (1LL << _width);
            _maxVal = (1LL << (_width - 1)) - 1;
        } else {
            _maxVal = 0x0001 << _width;
        }
    }

    int64_t BitFieldValue::value() const {
        if (_isSigned) {
		// Handle signed values
        	int64_t val = (_b & _mask) >> _offset;
        	if ((val & (1LL << (_width - 1))) != 0) {
        		val -= (1LL << _width);
        	}
        	return val;
        } else {
		// Handle unsigned values
        	return (_b & _mask) >> _offset;
        }
    }

    BitFieldValue& BitFieldValue::operator=(int64_t in) {
        // Error handle
	if (in < _minVal || in > _maxVal) {
            throw std::runtime_error("BitFieldValue: value out of range");
        }
        _b &= ~_mask;
        _b |= ((in << _offset) & _mask);
        return *this;
    }

    BitField64::BitField64(const std::string& initString) : _value(0), _joined(0) {
        init(initString);
    }

    BitField64::~BitField64() {
        for (auto field : _fields) {
            delete field;
        }
    }

    BitFieldValue& BitField64::operator[](size_t theIndex) {
        return *_fields.at(theIndex);
    }

    const BitFieldValue& BitField64::operator[](size_t theIndex) const {
        return *_fields.at(theIndex);
    }

    void BitField64::setFieldValue(const std::string& fieldName, int64_t value) {
        auto it = _map.find(fieldName);
        if (it == _map.end()) {
            throw std::runtime_error("BitField64: unknown field name");
        }
 
        unsigned int fieldIndex = it->second;
        (*this)[fieldIndex] = value;
    }

    size_t BitField64::index(const std::string& name) const {
        auto it = _map.find(name);
        if (it != _map.end()) {
            return it->second;
        } else {
            throw std::runtime_error("BitField64: unknown field name");
        }
    }

    BitFieldValue& BitField64::operator[](const std::string& name) {
        return *_fields[index(name)];
    }

    const BitFieldValue& BitField64::operator[](const std::string& name) const {
        return *_fields[index(name)];
    }

    unsigned BitField64::highestBit() const {
        unsigned hb = 0;
        for (const auto& field : _fields) {
            if (hb < (field->offset() + field->width())) {
                hb = field->offset() + field->width();
            }
        }
        return hb;
    }

    std::string BitField64::valueString() const {
        std::stringstream os;
        for (size_t i = 0; i < _fields.size(); ++i) {
            if (i != 0) os << ",";
            os << _fields[i]->name() << ":" << _fields[i]->value();
        }
        return os.str();
    }

    std::string BitField64::fieldDescription() const {
        std::stringstream os;
        for (size_t i = 0; i < _fields.size(); ++i) {
            if (i != 0) os << ",";
            os << _fields[i]->name() << ":" << _fields[i]->offset() << ":";
            if (_fields[i]->isSigned()) os << "-";
            os << _fields[i]->width();
        }
        return os.str();
    }

    void BitField64::addField(const std::string& name, unsigned offset, int width) {
        auto* bfv = new BitFieldValue(_value, name, offset, width);
        _fields.push_back(bfv);
        _map[name] = _fields.size() - 1;

        if (_joined & bfv->mask()) {
            throw std::runtime_error("BitField64::addField: bits already used");
        }

        _joined |= _fields.back()->mask();
    }

    void BitField64::init(const std::string& initString) {
	// Setup new objects
        unsigned offset = 0;
        std::vector<std::string> fieldDescriptors;
        std::istringstream ss(initString);
        std::string token;

	// Loop through the string and separte the values
        while (std::getline(ss, token, ',')) {
            fieldDescriptors.push_back(token);
        }

	// Set up each value separately
        for (const auto& desc : fieldDescriptors) {
            std::vector<std::string> subfields;
            std::istringstream ss2(desc);
            std::string subtoken;
	    // Look for subfields
            while (std::getline(ss2, subtoken, ':')) {
                subfields.push_back(subtoken);
            }

            std::string name;
            int width;
            unsigned thisOffset;
	    
	    // Handle subfields
            switch (subfields.size()) {
                case 2:
                    name = subfields[0];
                    width = std::stoi(subfields[1]);
                    thisOffset = offset;
                    offset += abs(width);
                    break;
                case 3:
                    name = subfields[0];
                    thisOffset = std::stoi(subfields[1]);
                    width = std::stoi(subfields[2]);
                    offset = thisOffset + abs(width);
                    break;
                default:
                    throw std::runtime_error("BitField64: invalid number of subfields");
            }

            addField(name, thisOffset, width);
        }
    }

    std::ostream& operator<<(std::ostream& os, const BitField64& b) {
        os << "bitfield: 0x" << std::hex << b.getValue() << std::dec << std::endl;
        for (const auto& it : b.getIndexMap()) {
            os << "  " << it.first << " [" << b[it.second].offset() << ":";
            if (b[it.second].isSigned()) os << "-";
            os << b[it.second].width() << "] : " << b[it.second].value() << std::endl;
        }
        return os;
    }
