#ifndef BITFIELD64_H
#define BITFIELD64_H

// Standard
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <cmath>
#include <algorithm>

/**
 * @class BitFieldValue
 * @brief Represents a value within a bitfield.
 * This is adapted from LCIO's BitFieldValue to work
 * for edm4hep. It is almost exactly the same.
 *
 * @author Samuel Ferraro
 */
class BitFieldValue {
    public:
	/**
         * @brief Constructor for BitFieldValue.
         * @param bitfield Reference to the containing bitfield.
         * @param name Name of the field.
         * @param offset Offset of the field within the bitfield.
         * @param signedWidth Width of the field. If negative, the field is signed.
         */	
        BitFieldValue(int64_t& bitfield, const std::string& name, unsigned offset, int signedWidth);
        /**
         * @brief Gets the value of the field.
         * @return The value of the field.
         */
	int64_t value() const;
        /**
         * @brief Sets the value of the field.
         * @param in The value to set.
         * @return Reference to the current object.
         */
	BitFieldValue& operator=(int64_t in);
        /**
         * @brief Implicit conversion to int64_t.
         * @return The value of the field.
         */
	operator int64_t() const { return value(); }
        /**
         * @brief Gets the name of the field.
         * @return The name of the field.
         */
	const std::string& name() const { return _name; }
        /**
         * @brief Gets the offset of the field.
         * @return The offset of the field.
         */
	unsigned offset() const { return _offset; }
        /**
         * @brief Gets the width of the field.
         * @return The width of the field.
         */
	unsigned width() const { return _width; }
        /**
         * @brief Checks if the field is signed.
         * @return True if the field is signed, false otherwise.
         */
	bool isSigned() const { return _isSigned; }
        /**
         * @brief Gets the bitmask of the field.
         * @return The bitmask of the field.
         */
	uint64_t mask() const { return _mask; }

    private:
	int64_t& _b;            ///< Reference to the containing bitfield.
        uint64_t _mask;         ///< Bitmask for the field.
        std::string _name;      ///< Name of the field.
        unsigned _offset;       ///< Offset of the field within the bitfield.
        unsigned _width;        ///< Width of the field.
        int _minVal;            ///< Minimum value for the field (if signed).
        int _maxVal;            ///< Maximum value for the field (if signed).
        bool _isSigned;         ///< Flag indicating if the field is signed.



    };

/**
 * @class BitField64
 * @brief Represents a 64-bit bitfield with multiple fields.
 * This is adapted from LCIO's BitField64 to work with
 * edm4hep. It is almost exactly the same.
 *
 * @author Samuel Ferraro
 */
class BitField64 {
    public:
        using IndexMap = std::map<std::string, unsigned int>;
	
	/**
         * @brief Constructor for BitField64.
         * @param initString Initialization string describing the fields.
         */
        BitField64(const std::string& initString);
        ~BitField64();
        /**
         * @brief Gets the value of the bitfield.
         * @return The value of the bitfield.
         */
	int64_t getValue() const { return _value; }
        /**
         * @brief Sets the value of the bitfield.
         * @param value The value to set.
         */
	void setValue(int64_t value) { _value = value; }
        /**
 	 * @brief Sets the value of one of the components of the bitfield.
 	 * @param fieldName The string associated with the field.
 	 * @param value The value to fill the field with
 	 */
	void setFieldValue(const std::string& fieldName, int64_t value);
	/**
         * @brief Resets the bitfield to zero.
         */
	void reset() { _value = 0; }
        /**
         * @brief Accesses a field by index.
         * @param theIndex The index of the field.
         * @return Reference to the BitFieldValue.
         */
	BitFieldValue& operator[](size_t theIndex);
        /**
         * @brief Accesses a field by index (const version).
         * @param theIndex The index of the field.
         * @return Const reference to the BitFieldValue.
         */
	const BitFieldValue& operator[](size_t theIndex) const;
        /**
         * @brief Gets the highest bit set in the bitfield.
         * @return The highest bit set.
         */
	unsigned highestBit() const;
        /**
         * @brief Gets the number of fields in the bitfield.
         * @return The number of fields.
         */
	size_t size() const { return _fields.size(); }
        /**
         * @brief Gets the index of a field by name.
         * @param name The name of the field.
         * @return The index of the field.
         */
	size_t index(const std::string& name) const;
        /**
         * @brief Accesses a field by name.
         * @param name The name of the field.
         * @return Reference to the BitFieldValue.
         */
	BitFieldValue& operator[](const std::string& name);
        /**
         * @brief Accesses a field by name (const version).
         * @param name The name of the field.
         * @return Const reference to the BitFieldValue.
         */
	const BitFieldValue& operator[](const std::string& name) const;
        /**
         * @brief Gets the low word (lower 32 bits) of the bitfield.
         * @return The low word.
         */
	unsigned lowWord() const { return static_cast<unsigned>(_value & 0xffffFFFF); }
        /**
         * @brief Gets the high word (upper 32 bits) of the bitfield.
         * @return The high word.
         */
	unsigned highWord() const { return static_cast<unsigned>(_value >> 32); }
        /**
         * @brief Gets a description of the fields in the bitfield.
         * @return A string describing the fields.
         */
	std::string fieldDescription() const;
        /**
         * @brief Gets a string representation of the bitfield values.
         * @return A string representation of the bitfield values.
         */
	std::string valueString() const;
        /**
	 * @brief Gets the index map for looping
	 * @return An Index Map
	 */
	const IndexMap& getIndexMap() const { return _map; }

    protected:
	/**
         * @brief Adds a field to the bitfield.
         * @param name The name of the field.
         * @param offset The offset of the field.
         * @param width The width of the field.
         */
        void addField(const std::string& name, unsigned offset, int width);
        /**
         * @brief Initializes the bitfield using an initialization string.
         * @param initString The initialization string.
         */
	void init(const std::string& initString);

	std::vector<BitFieldValue*> _fields;  ///< Vector of BitFieldValue pointers.
        int64_t _value;                        ///< The value of the bitfield.
        IndexMap _map;                        ///< Map of field names to indices.
        uint64_t _joined;                      ///< Bitmask representing the used bits.

    };

    /**
     * @brief Output stream operator for BitField64.
     * @param os The output stream.
     * @param b The BitField64 object.
     * @return The output stream.
     */
    std::ostream& operator<<(std::ostream& os, const BitField64& b);
#endif //BITFIELD64_H
